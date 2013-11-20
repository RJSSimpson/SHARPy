'''@package Linear dynamic simulations of coupled systems.

Created on 10 Sep 2013
@author: rjs10
'''

from PyMPC.ssdiscrete import StateSpace
from scipy.signal import dlsim
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.io import loadmat
from scipy.interpolate import interp1d
import SharPySettings as Settings
from collections import OrderedDict
from PyAero.UVLM.Utils.DerivedTypesAero import Gust
from PyMPC import MPC

# Flag for some interactive routines
Interactive = True

def Solve_Py(*args, **kwords):
    """@brief time domain solution of linear aeroelastic system.
    args:
        * 4 (matPath/(A,B,C,D,dt), U, t, x0)
        
    @param matPath string containing absolute path of .mat file with disSys
            struct.
    @param A state space system matrix.
    @param B state space system matrix.
    @param C state space system matrix.
    @param D state space system matrix.
    @param dt Timestep.
    @param U u(t), the input sequence (None for free response).
    @param t timesteps at which input is defined (None to default).
    @param x0 The initial state size(x,1) (None: zero by default).
    @param writeDict OrderedDict of 'name':output index to write in loop, Keyword.
    @param plotDict OrderedDict of 'name':output index to plot in loop, Keyword.
    @returns tout Time values for output.
    @returns yout y(t), the system response.
    @returns xout x(t), the state history.
    """
    
    if len(args) == 4:
        # Initialize discrete time system
        if isinstance(args[0], str):
            disSys = StateSpace(args[0])
        elif isinstance(args[0], tuple):
            disSys = StateSpace(args[0][0],
                                args[0][1],
                                args[0][2],
                                args[0][3],
                                disSys.A)
        elif isinstance(args[0],StateSpace):
            disSys = args[0]
        else:
            raise TypeError("First argument (of 4) not recognised.")
        
        # check input sequence
        try:
            disSys._Ts.shape[1] #check array has 2 dimensions
            arg2d = args[1] 
        except:
            arg2d = np.atleast_2d(args[1]) #make 2d for comparison with nU
            arg2d = arg2d.T
            
        # assign time to local var
        t = args[2]
        
        if (args[1] == None or 'mpcCont' in kwords) and t == None:
            if Interactive == True:
                userNumber = input('Enter number of time steps (dt = %f):'\
                                   % disSys._Ts)
                userNumber = int(userNumber)
                assert (userNumber > 0 and userNumber < 1000000),\
                       IOError("Number of steps must be > 0 and < 10^6.")
            else:
                userNumber = np.ceil(1.0/disSys.dt)
            
            U = np.zeros((userNumber, disSys.nU))
        elif disSys.nU != arg2d.shape[1]:
            raise ValueError("Wrong number of inputs in input sequence")
        else:
            U = arg2d
            
        # Run simulation
        if 'writeDict' in kwords or 'plotDict' in kwords \
        or 'mpcCont' in kwords or 'gust' in  kwords:
            # Determine whether to write and/or plot
            if 'writeDict' in kwords and Settings.WriteOut == True:
                write = True
            
            if 'plotDict' in kwords and Settings.PlotOut == True:
                #plot = True
                raise NotImplementedError()
            
            # Check function inputs for discrete time system plotting/writing
            if t is None:
                out_samples = U.shape[0]
                stoptime = (out_samples - 1) * disSys._Ts
            else:
                stoptime = t[-1]
                out_samples = int(np.floor(stoptime / disSys._Ts)) + 1
        
            # Pre-build output arrays
            xout = np.zeros((out_samples, disSys.A.shape[0]))
            yout = np.zeros((out_samples, disSys.C.shape[0]))
            tout = np.linspace(0.0, stoptime, num=out_samples)
        
            # Check initial condition
            if args[3] is None:
                xout[0,:] = np.zeros((disSys.A.shape[1],))
            else:
                xout[0,:] = np.asarray(args[3])
        
            # Pre-interpolate inputs into the desired time steps
            if t is None and U[0,0] is not None:
                u_dt = U
            elif U[0,0] is None and 'mpcCont' in kwords:
                u_dt = np.zeros((t.shape[0],disSys.nU))
            else:
                if len(U.shape) == 1:
                    U = U[:, np.newaxis]
        
                u_dt_interp = interp1d(t, U.transpose(),
                                       copy=False,
                                       bounds_error=True)
                u_dt = u_dt_interp(tout).transpose()
                
            if write == True:
                # Write output file header
                outputIndices = list(writeDict.values())
                ofile = Settings.OutputDir + \
                        Settings.OutputFileRoot + \
                        '_SOL302_out.dat'
                fp = open(ofile,'w')
                fp.write("{:<14}".format("Time"))
                for output in writeDict.keys():
                    fp.write("{:<14}".format(output))
                fp.write("\n")
                fp.flush()
            # END if write
                
            # Simulate the system
            for i in range(0, out_samples - 1):
                
                # get optimal control action
                if 'mpcCont' in kwords:             
                    u_dt[i,:kwords['mpcCont'].mpcU] = kwords['mpcCont'].getUopt(xout[i],True)
                    
                # get gust velocity at current time step
                if 'gust' in kwords:
                    raise NotImplementedError()
                
                xout[i+1,:] = np.dot(disSys.A, xout[i,:]) + np.dot(disSys.B, u_dt[i,:])
                yout[i,:] = np.dot(disSys.C, xout[i,:]) + np.dot(disSys.D, u_dt[i,:])
                
                if write == True:
                    fp.write("{:<14,e}".format(tout[i]))
                    for j in yout[i,outputIndices]:
                        fp.write("{:<14,e}".format(j))
                    fp.write("\n")
                    fp.flush()
        
            # Last point
            yout[out_samples-1,:] = np.dot(disSys.C, xout[out_samples-1,:]) + \
                                    np.dot(disSys.D, u_dt[out_samples-1,:])
            # Write final output and close
            fp.write("{:<14,e}".format(tout[out_samples-1]))
            for j in yout[out_samples-1,outputIndices]:
                fp.write("{:<14,e}".format(j))      
            fp.close()
                                    
            return tout, yout, xout
        else:
            # Run discrete time sim using scipy solver
            tout, yout, xout = dlsim((disSys.A,disSys.B,
                                      disSys.C,disSys.D,
                                      disSys._Ts),
                                      U,
                                      t = args[2],
                                      x0 = args[3])
        return tout, yout, xout
    else:
        raise ValueError("Needs 4 positional arguments")
    
def plotOutputs(tout, yout, xout, youtJ = None, xoutJ = None):
    """@brief plot outputs from linear system.
    @param tout time series data.
    @param yout output series data.
    @param xout state series data.
    @param youtJ indices of outputs of interest (None = plot all).
    @param xoutJ indices of states of interest (None = plot none).
    """
    if youtJ == None:
        nY = yout.shape[1]
        # fin number of columns in figure
        if nY <= 3:
            nC = nY
        elif nY > 3:
            nC = 3
        for j in range(nY):
            plt.subplot(np.ceil(nY/3.0),nC,j+1) # plot each output
            plt.ylabel(r'$y_{%d}$' % j)
            plt.xlabel(r'$t [secs]$')
            plt.plot(tout[:],yout[:,j],'ko-')
    elif isinstance(youtJ,tuple):
        nY = len(youtJ)
        # fin number of columns in figure
        if nY <= 3:
            nC = nY
        elif nY > 3:
            nC = 3 
        # initalise plot counter
        plotCount = 0
        for j in youtJ:
            plt.subplot(np.ceil(nY/3.0),nC,plotCount+1) # plot each output
            plt.ylabel(r'$y_{%d}$' % j)
            plt.xlabel(r'$t [secs]$')
            plt.plot(tout[:],yout[:,j],'ko-')
            plotCount += 1
            
    else:
        raise TypeError("youtJ must be None or a tuple of indices")
    plt.show()
    if xoutJ != None:
        raise NotImplementedError("plotting of states not implemented yet.")
    
def removeGustInputs(disSys):
    """@brief Return state-space system with gust inputs removed.
    @param disSys Discrete-time state-space object."""
    # Get number of panels to determine number of gust inputs
    assert os.path.isfile(disSys._matPath + '.mat'), \
                   IOError("File doesn't exist")
    myDict = loadmat(disSys._matPath, variable_names = ['Mwing','Nwing'])
    nUg = 3*myDict['Mwing'][0,0]*myDict['Nwing'][0,0] # Number of gust inputs
    # Check to see if nU is greater than nUg
    if disSys.nU <=  nUg:
        raise ValueError("System inputs is less than number of gust inputs.")
    else:
        nUc = disSys.nU - nUg
        return StateSpace(disSys.A,
                          disSys.B[:,:nUc],
                          disSys.C,
                          np.zeros([disSys.nY,nUc]),
                          disSys.dt)
    
if __name__ == '__main__':
    
    matPath = '/home/rjs10/git/SHARP/output/Goland/TorsionBending/Q140_N8_Py'
    # Test remove gust inputs
    disSys = StateSpace(matPath)
    #disSys = removeGustInputs(disSys)
    # Define 1 degree sinusoidal control input
    betaHat = 1.0*np.pi/180.0
    omega = 30.0 # rads-1
    t = np.arange(0,0.5,disSys.dt)
    U = betaHat*np.sin(omega*t)
    # Define names and indices of outputs.
    writeDict = OrderedDict()
    writeDict['kappa_x_root'] = 0
    writeDict['kappa_y_root'] = 1
    # Initialize gusts.
    gust = Gust(uMag = 0.01*140.0,
                l = 10.0,
                r = 0.0)
    # Initialize MPC controller.
    mpcCont = MPC.MPC('golandControl','/home/rjs10/git/SHARPy/src/PyMPC/systems/')
    # Initial conditions
    x0 = np.zeros((disSys.A.shape[0]))
    x0[0] = 1
    # Run solver
    tout, yout, xout = Solve_Py(disSys, None, t, x0,
                                writeDict = writeDict,
                                mpcCont = mpcCont)
    plotOutputs(tout, yout, xout, youtJ = None)
    