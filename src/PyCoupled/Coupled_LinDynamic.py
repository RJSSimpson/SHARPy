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
        * 4 (matPath, U, t, x0)
        * 5(+) NotImplemented (TODO: write a linearization routine.)
        
    @param matPath string containing absolute path of .mat file with disSys struct.
    @param U u(t), the input sequence. can be zeros to simulate free response.
    @param t timesteps at which input is defined (None to default).
    @param x0 The initial condition on the state vector(None: zero by default).
    @returns tout Time values for output.
    @returns yout y(t), the system response.
    @returns xout x(t), the state history.
    """
    
    if len(args) == 4:
        # Initialize discrete time system
        disSys = StateSpace(args[0])
        # check input sequence
        if disSys.nU() != args[1].shape[1]:
            raise ValueError("Wrong number of inputs in input sequence")
        
        # Run discrete time sim
        dlsim((disSys.A,disSys.B,disSys.C,disSys.D,disSys._Ts),
              U,
              t = args[2],
              x0 = args[3])
        
    elif len(args) >= 5:
        raise NotImplementedError("Linearization not available in SHARPy yet.")
    
    
if __name__ == '__main__':
    
    matPath = '/home/rjs10/Documents/MATLAB/MPC/HALE/HALE_sigma1000_Py'
    U = np.zeros((2,4)) # HALE example, no thrust
    
    Solve_Py(matPath,U,None,None)