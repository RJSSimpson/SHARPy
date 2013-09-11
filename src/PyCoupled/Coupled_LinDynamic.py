'''@package Linear dynamic simulations of coupled systems.

Created on 10 Sep 2013
@author: rjs10
'''

from PyMPC.ssdiscrete import StateSpace
from scipy.signal import dlsim
import numpy as np

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
    @param U u(t), the input sequence. can be zeros to simulate free response.
    @param t timesteps at which input is defined (None to default).
    @param x0 The initial condition on the state vector(None: zero by default).
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
        # check input sequence
        if disSys.nU() != args[1].shape[1]:
            raise ValueError("Wrong number of inputs in input sequence")
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
    
    
if __name__ == '__main__':
    
    matPath = '/home/rjs10/Documents/MATLAB/MPC/HALE/HALE_sigma1000_Py'
    U = np.zeros((1,4)) # HALE example, no thrust
    
    Solve_Py(matPath,U,None,None)
    
    # Define system
    A=np.array([[1]])
    B=np.array([[1,1,1]])
    C=np.array([[1]])
    D=np.array([[0,0,0]])
    dt = 1.0
    
    # Define input
    U=np.array([[0,0,0]]) # Input tensor with 1 timestep.
    
    dlsim((A,B,C,D,dt),U)
    
    Solve_Py((A,B,C,D,dt),U,None,None)