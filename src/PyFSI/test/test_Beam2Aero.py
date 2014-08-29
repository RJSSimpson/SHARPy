'''PyFSI.test.test_Beam2Aero
@brief      Unit tests for mapping from beam to aero.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       17/07/2014
@pre        None
@warning    None
'''
import unittest
import SharPySettings as Settings
import XbeamLib as xbl
import numpy as np
from scipy.linalg import logm, expm
from PyFSI.Beam2Aero_Lin import zetaDotSubMat

TestDir = Settings.SharPyProjectDir + 'output/tests/PyFSI/Beam2Aero/'

class Test_NonlinearGridVelocity(unittest.TestCase):
    """@brief Compare velocities from different formulations."""

    def setUp(self):
        """@brief Set SharPy output directory and file root."""
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = ''


    def tearDown(self):
        pass


    def test_CompareVels(self):
        """Calculate using nonlinear formulations and compare.
        @param a0 displacement amplitude scale factor.
        @param a1 rotation amplitude scale factor.
        @param rot rotation matrix from G to a frame.
        @param rotDot time derivative of rot.
        """
        
        # parameters
        off0 = 1.0
        off1 = 1.0*np.pi/2
        a0=1.0
        a1=1.0*np.pi/2
        omega=2*np.pi
        
        # time histories
        t = np.linspace(0.0, 2.0, 200)
        v_a = np.zeros((3,len(t)))
        omega_a=np.zeros((3,len(t)))
        r=np.zeros((3,len(t)))
        rDot=np.zeros((3,len(t)))
        psi=np.zeros((3,len(t)))
        psiDot=np.zeros((3,len(t)))
        xi=np.zeros((3,len(t)))
        xiDot=np.zeros((3,len(t)))
        psi_G2a=np.zeros((3,len(t)))
        velNonlin=np.zeros((3,len(t)))
        velLin=np.zeros((3,len(t)))
        
        # initial states
        v_a[:,0]=off0+a0*np.sin(omega*t[0])
        omega_a[:,0]=off1+a1*np.sin(omega*t[0])
        r[:,0]=off0+a0*np.sin(omega*t[0])
        rDot[:,0]=omega*a0*np.cos(omega*t[0])
        psi[:,0]=off1+a1*np.sin(omega*t[0])
        psiDot[:,0]=omega*a1*np.cos(omega*t[0])
        xi[:,0]=off0+a0*np.sin(omega*t[0])
        xiDot[:,0]=omega*a0*np.cos(omega*t[0])
        
        # rotation matrix for integrating omega_a
        rot = np.zeros((3,3))
        rotDot = np.zeros((3,3))
        
        # column matrix of structural DoFs
        nBeam=1 #one beam node
        nSection=1 #one aero node per section
        nStates = 2*6*nBeam + 9 + 2*3*nSection
        x = np.zeros((nStates))
        
        # generate linear function
        Xi = zetaDotSubMat(r[:,0],psi[:,0],rDot[:,0],psiDot[:,0],v_a[:,0],
                           omega_a[:,0],psi_G2a[:,0],xi[:,0],xiDot[:,0],0,0,1,1)
        
        # initialise time histories
        for j in np.arange(len(t)):
            # create latest states
            v_a[:,j]=off0+a0*np.sin(omega*t[j])
            omega_a[:,j]=off1+a1*np.sin(omega*t[j])
            r[:,j]=off0+a0*np.sin(omega*t[j])
            rDot[:,j]=omega*a0*np.cos(omega*t[j])
            psi[:,j]=off1+a1*np.sin(omega*t[j])
            psiDot[:,j]=omega*a1*np.cos(omega*t[j])
            xi[:,j]=off0+a0*np.sin(omega*t[j])
            xiDot[:,j]=omega*a0*np.cos(omega*t[j])
            
            # integrate omega_a to get psi_G2a
            rot = expm(xbl.Skew(psi_G2a[:,j]))
            rotDot = np.dot(rot,xbl.Skew(omega_a[:,j]))
            if j > 0:
                rot = rot + (t[j]-t[j-1])*rotDot
                psi_G2a[:,j] = xbl.VectofSkew(np.real(logm(rot)))
            # end if j > 0
            
            # calculate nonlinear velocity
            velNonlin[:,j]=NonlinearGridVels(v_a[:,j], omega_a[:,j],
                                             psi_G2a[:,j], r[:,j], rDot[:,j], 
                                             psi[:,j], psiDot[:,j], xi[:,j],
                                             xiDot[:,j])
            
            # calculate linear velocity
            # assign to column matrix
            x[:]=np.concatenate((r[:,j],psi[:,j],rDot[:,j],psiDot[:,j],v_a[:,j],
                                 omega_a[:,j],psi_G2a[:,j],xi[:,j],xiDot[:,j]))
            
            velLin[:,j]=np.dot(Xi,x)
            
            print('x = ', x)
            print('Lin - Nonlin = ',velLin[:,j]-velNonlin[:,j]) #TODO: this is giving epsilon size errors only!?
        # end for j in len(t)
        pass

def NonlinearGridVels(v_a,omega_a,psi_G2a,r,rDot,psi,psiDot,xi,xiDot):
    """@brief Calculate grid velocities based on pre-linearized (nonlinear)
    formulation.
    @param r Position vector on beam, defined in the a-frame.
    @param psi CRV of orientation of beam cross-section.
    @param rDot Velocity of beam at r, defined in the a-frame.
    @param psiDot Time rate of change of psi.
    @param v_a Velocity of a-frame, defined in the a-frame.
    @param omega_a Angular vel of a-frame, defined in a-frame.
    @param psi_G2a CRV of orientation of a-frame relative to earth.
    @param xi Cross-sectional coordinate, defined in B-frame.
    @param xiDot Time rate of change of xi, defined in B-frame.
    @returns gridVel velocity of point on aero surface in the inertial G-frame.
    @details Terms in the formulation below correspond to those in
    PyFSI.Beam2Aero_Lin.zetaDotSubMat.
    """
    gridVel = np.zeros((3))
    # transformation matrices
    cGa = xbl.Psi2TransMat(psi_G2a)
    skewOmAa = xbl.Skew(omega_a)
    cAb = xbl.Psi2TransMat(psi)
    tangAb = xbl.Tangential(psi)
    #term 1:
    gridVel[:] += np.dot(cGa,v_a)
    #term 2:
    gridVel[:] += np.dot(cGa,
                  np.dot(skewOmAa,
                         r))
    # term 3:
    gridVel[:] += np.dot(cGa,
                  np.dot(skewOmAa,
                  np.dot(cAb,
                         xi)))
    #terms 4:
    gridVel[:] += np.dot(cGa,rDot)
    #term 5:
    gridVel[:] += np.dot(cGa,
                  np.dot(cAb,
                         xiDot))
    #term 6:
    gridVel[:] += np.dot(cGa,
                  np.dot(cAb,
                  np.dot(xbl.Skew(np.dot(tangAb,psiDot)),
                  xi)))
    return gridVel
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()