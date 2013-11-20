'''
Created on 1 Nov 2013

@author: rjs10
'''

import sys
import muaompc
import numpy as np
import os
from scipy.io import loadmat

class MPC:
    """@brief MPC functions."""
    
    def __init__(self, moduleName, moduleDir, LQR = False):
        """@brief Control law formulation and code generation.
        
        @param moduleName Name of module to import control specification from.
        @param moduleDir Directory containing module.
        @param LQR Flag to use simple constant gain LQR control instead.
        """
        
        self.moduleName = moduleName
        self.moduleDir = moduleDir
        
        # Load the module.
        sys.path.append(moduleDir)
        mod = __import__(moduleName)
        self.Ad = mod.Ad
        self.Bd = mod.Bd
        
        # Save attributes.
        self.T = mod.T
        self.nFull = self.T.shape[1] # Full-order linear system states.
        self.mpcU = self.Bd.shape[1] # Number of MPCed inputs
        self.nAux = mod.nAux
        self.matPath = mod.matPath
        
        # Do LQR?
        self.LQR = LQR
        if self.LQR == True:
            self.u_lb = mod.u_lb
            self.u_ub = mod.u_ub
            self.e_lb = mod.e_lb[-self.nAux:] # mixed constraints corresponding
            self.e_ub = mod.e_ub[-self.nAux:] # to inputs rates.
            # get gain matrix
            matDir = self.matPath.rsplit('/',1)
            kPath = matDir[0] + "/Q140_FullSystem_GustInputs_L10_lqrCont_lqrKforPy"
            Dict = loadmat(kPath,variable_names = ['K'])
            K = np.zeros((self.mpcU,self.Ad.shape[0]))
            K[:,:-self.nAux] = Dict['K'].copy('C')
            self.K = K
            
            # data to saturate input rates
            self.u_sav = 0.0
            
        else:
            # Generate MPC control law and files.
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
                self.ctl.x_ref[self.mpc.size.states*i:\
                               self.mpc.size.states*(i+1)]=\
                               self.setPoint
    
    def getUopt(self,xFull,xRedFlag = False):
        """@brief get optimal control input for the next timestep.
        
        @param xFull Full-order state measurement.
        @param xRedFlag True if the state is already reduced.
        """
        if xRedFlag == False:
            xRed = np.dot(self.T,xFull) # Get reduced state.
        else:
            xRed = np.zeros((xFull.shape[0] + self.nAux))
            xRed[:xFull.shape[0]] = xFull.copy('C')
            
        if self.LQR == True:
            u_opt = -np.dot(self.K,xRed)
            # apply saturation
#             if u_opt[0] < self.u_lb[0][0]:
#                 u_opt[0] = self.u_lb[0][0]
#             elif u_opt[0] > self.u_ub[0][0]:
#                 u_opt[0] = self.u_ub[0][0]
                
            # apply rate saturation
#             d_u_opt = u_opt - self.u_sav
#             if d_u_opt < self.e_lb[0][0]:
#                 u_opt[0] = self.u_sav + self.e_lb[0][0] 
#             elif d_u_opt > self.e_ub[0][0]:
#                 u_opt[0] = self.u_sav + self.e_ub[0][0]
            
            # save and return as 2d array
            self.u_sav = u_opt
            u_2darray = np.zeros((1,1))
            u_2darray[0,0] = u_opt[0]
            return(u_2darray)
        else:
            xDist = xRed - self.x_opt_kp1 # Calculate disturbance based on previous estimate.
            
            # Apply disturbance along prediction horizon.
    #         for i in range(self.mpc.N):
    #             self.ctl.x_ref[self.mpc.size.states*i:self.mpc.size.states*(i+1),0] =\
    #                                                         (self.setPoint.flatten() - xDist)
                                                            
            # Solve optimization problem.
            self.ctl.solve_problem(xRed)
            
            # Predict expected output at next time step.
            self.x_opt_kp1 = np.dot(self.Ad,xRed) + np.dot(self.Bd,self.ctl.u_opt[:self.mpc.size.inputs].flatten())
            
            # Return control action.
            return(self.ctl.u_opt[:self.mpc.size.inputs])
    
    def genForMatlab(self,matDir=None):
        """@brief Generate files for using muao-mpc control in matlab.
        @param matDir The path to where the cfiles should be copied."""
        
        # Load the module.
        sys.path.append(self.moduleDir)
        
        # Generate Files
        mpcMat = muaompc.ltidt.setup_mpc_problem(self.moduleName)
        mpcMat.generate_c_files(matlab=True)
        print(self.matPath)
        if matDir == None:
            os.system("cp -R cmpc/ " + self.matPath + "_cmpc")
        else:
            os.system("cp -R cmpc/ " + matDir)

if __name__ == '__main__':
    myMpc = MPC('golandControl','/home/rjs10/git/SHARPy/src/PyMPC/systems/', LQR = True)
    #myMpc.genForMatlab() # Generate file for matlab