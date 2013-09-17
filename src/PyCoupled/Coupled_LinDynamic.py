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
                                args[0][4])
        elif isinstance(args[0],StateSpace):
            disSys = args[0]
        else:
            raise TypeError("First argument (of 4) not recognised.")
        # check input sequence
        try:
            args[1].shape[1] #check array has 2 dimensions
            arg2d = args[1] 
        except:
            arg2d = np.atleast_2d(args[1]) #make 2d for comparison with nU
            arg2d = arg2d.T
        if args[1] == None:
            if Interactive == True:
                userNumber = input('Enter number of time steps (dt = %f):'\
                                   % disSys._Ts)
                userNumber = int(userNumber)
                assert (userNumber > 0 and userNumber < 1000000),\
                       IOError("Number of steps must be > 0 and < 10^6.")
            else:
                userNumber = 100
            
            U = np.zeros((userNumber, disSys.nU))
        elif disSys.nU != arg2d.shape[1]:
            raise ValueError("Wrong number of inputs in input sequence")
        else:
            U = arg2d
        # Run discrete time sim
        tout, yout, xout = dlsim((disSys.A,disSys.B,
                                  disSys.C,disSys.D,
                                  disSys._Ts),
                                  U,
                                  t = args[2],
                                  x0 = args[3])
        return tout, yout, xout
    else:
        raise ValueError("Needs 4 input arguments")
    
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
        raise TypeError("youJ must be None or a tuple of indices")
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
    
    matPath = '/home/rjs10/git/SHARP/output/Goland/Goland_ss_Q140_M16N28_wake27_ROM_Py'
    # Test remove gust inputs
    disSys = StateSpace(matPath)
    disSysBc = removeGustInputs(disSys)
    # Define 1 degree sinusoidal control input
    betaHat = 1.0*np.pi/180.0
    omega = 30.0 # rads-1
    t = np.arange(0,0.5,disSysBc.dt)
    U = betaHat*np.sin(omega*t)
    tout, yout, xout = Solve_Py(disSysBc,U,None,None)
    plotOutputs(tout, yout, xout, youtJ = None)
    