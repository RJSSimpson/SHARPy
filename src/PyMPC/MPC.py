'''
Created on 1 Nov 2013

@author: rjs10
'''

import sys
import muaompc
import numpy as np

class MPC:
    """@brief MPC functions."""
    
    def __init__(self, moduleName, moduleDir):
        """@brief Control law formulation and code generation.
        
        @param modulePath Full path of module to import control specification
                from.
        """
        
        # Load the module.
        sys.path.append(moduleDir)
        mod = __import__(moduleName)
        self.Ad = mod.Ad
        self.Bd = mod.Bd
        
        # Save attributes.
        self.T = mod.T
        self.nFull = self.T.shape[1] # Full-order linear system states.
        
        # Generate control law and files.
        self.mpc = muaompc.ltidt.setup_mpc_problem(moduleName)
        self.mpc.generate_c_files()
    
        # Configure optimizer.
        self.ctl = self.mpc.ctl
        self.ctl.conf.in_iter = 24 # Internal iterations in augmented Lagr method
        self.ctl.conf.ex_iter  =2  # External iteration in augmented Lagr method
        self.ctl.conf.warmstart = True
        
        # initialise prediction from previous timestep as zero
        self.x_opt_kp1 = np.zeros(self.Ad.shape[0])
        
        # Initialize undisturbed reference trajectory
        self.setPoint = np.zeros((self.Ad.shape[0],1))    
        for i in range(self.mpc.N):
            self.ctl.x_ref[self.mpc.size.states*i:self.mpc.size.states*(i+1)] =\
                                                                   self.setPoint
    
    def getUopt(self,xFull):
        """@brief get optimal control input for the next timestep.
        
        @param xFull Full-order state measurement.
        """
        xRed = np.dot(self.T,xFull) # Get reduced state.
        xDist = xRed - self.x_opt_kp1 # Calculate disturbance basedon previous estimate.
        
        # Apply disturbance along prediction horizon.
#         for i in range(self.mpc.N):
#             self.ctl.x_ref[self.mpc.size.states*i:self.mpc.size.states*(i+1),0] =\
#                                                         (self.setPoint.flatten() - xDist)#*np.exp(-1.0*(i/(self.mpc.N-1)))
                                                        
        # Solve optimization problem.
        self.ctl.solve_problem(xRed)
        
        # Predict expected output at next time step.
        self.x_opt_kp1 = np.dot(self.Ad,xRed) + np.dot(self.Bd,self.ctl.u_opt[:self.mpc.size.inputs].flatten())
        
        # Return control action.
        return(self.ctl.u_opt[:self.mpc.size.inputs])
        

if __name__ == '__main__':
    myMpc = MPC('golandControl','/home/rjs10/git/SHARPy/src/PyMPC/systems/')
    xFull = np.zeros((myMpc.nFull))
    print(myMpc.getUopt(xFull))