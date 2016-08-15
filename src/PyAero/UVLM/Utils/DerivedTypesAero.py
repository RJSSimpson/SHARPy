'''PyAero.UVLM.Utils.DerivedTypes
@brief      Classes containing data for UVLM solver.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       30/01/2013
@pre        None
@warning    None
'''

import ctypes as ct
import numpy as np
from PyAero.UVLM.Utils.UVLMLib import Colloc

class VMopts:
    """@brief options for UVLM
    @param M Chordwise panels.
    @param N Spanwise panels.
    @param ImageMethod Use image method across y-z plane if True.
    @param Mstar Number of wake panels.
    @param Steady Flag for specifying Steady(True) or Unsteady(False) solution.
    @param KJMeth If True use Joukowski force calculation method, if false use
           generalised Katz and Plotkin.
    @param NewAIC Calculate a new AIC matrix, required if the aerodynamic 
           surface has deformed significantly.
    @param DelTime Desired time step.
    @param Rollup Flag to specify wake roll-up due to induced velocities.
    @param NumCores Number of cores for parallel computation using openMP."""
    
    def __init__(self, M, N, ImageMethod = False, Mstar = 1, Steady = True,\
                 KJMeth = True, NewAIC = True, DelTime = 0.0, Rollup = False,\
                 NumCores = 1):
        self.M = ct.c_uint(M)
        self.N = ct.c_uint(N)
        self.ImageMethod = ct.c_bool(ImageMethod)
        self.Mstar = ct.c_int(Mstar)
        self.Steady = ct.c_bool(Steady)
        self.KJMeth = ct.c_bool(KJMeth)
        self.NewAIC = ct.c_bool(NewAIC)
        self.DelTime = ct.c_double(DelTime)
        self.Rollup = ct.c_bool(Rollup)
        self.NumCores = ct.c_uint(NumCores)
        

class VMinput:
    """@brief input data for UVLM.
    @param c chord.
    @param b span.
    @param area wing area.
    @param U_infty Free-stream flow magnitude.
    @param alpha free-stream direction (AoA).
    @param theta root twist angle.
    @param WakeLength length of wake in chords.
    @param ctrlSurf ControlSurf object.
    """
    
    def __init__(self, c, b, U_mag, alpha, theta,
                  WakeLength = 50.0,
                  ctrlSurf = None,
                  gust = None):
        self.c = c
        self.b = b
        self.area = c*b
        self.U_infty = np.zeros((3))
        self.U_infty[:] = [0.0, -U_mag*np.cos(alpha), U_mag*np.sin(alpha)]
        self.alpha = alpha
        self.theta = theta
        self.WakeLength = WakeLength
        self.ctrlSurf = ctrlSurf
        self.gust = gust
        
class ControlSurf:
    """@brief Control surface data and member functions.
    @param iMin index at which control surface starts in chordwise sense.
    @param iMax index at which control surface ends in chordwise sense.
    @param jMin index at which control surface starts in spanwise sense.
    @param jMax index at which control surface ends in spanwise sense.
    @param typeMotion string describing type of prescribed, or other, motion.
    @param betaBar Amplitude of control surface angle.
    @param omega Angular rate of control surface oscillations [rad/s].
    @details Control surface rates are calculated using backward difference
             in update and init routines.
    """
    
    def __init__(self, iMin, iMax, jMin, jMax,
                  typeMotion, betaBar,
                  omega = 0.0):
        """@brief Initialise control surface based on typeMotion."""
        
        self.iMin = iMin
        self.iMax = iMax
        self.jMin = jMin
        self.jMax = jMax
        self.typeMotion = typeMotion
        self.betaBar = betaBar
        self.omega = omega
        
        if typeMotion == 'asInput':
            self.beta = 0.0
            self.betaDot = 0.0
        elif typeMotion == 'sin':
            self.beta = 0.0
            self.betaDot = omega * betaBar
        elif typeMotion == 'cos':
            self.beta = betaBar
            self.betaDot = 0.0
        else:
            raise Exception('Control Surface typeMotion not recognised')
        
        # Define start time for calculation of derivative
        self.time = 0.0
        
        
    def update(self, time, beta = None):
        """@brief Calculate new beta and betaDot."""
        
        if (beta == None):
            if (self.typeMotion == 'sin'):
                self.beta = self.betaBar * np.sin(self.omega*time)
                self.betaDot = self.omega*self.betaBar*np.cos(self.omega*time)
                self.time = time # For homogeneous function of update routine 
            elif (self.typeMotion == 'cos'):
                self.beta = self.betaBar * np.cos(self.omega*time)
                self.betaDot = - self.omega*self.betaBar*np.sin(self.omega*time)
                self.time = time # For homogeneous function of update routine
            elif (self.typeMotion == 'asInput'):
                self.beta = 0.0
                self.betaDot = 0.0
                self.time = time
            else:
                raise Exception ('beta unspecified and typeMotion not ' + 
                                 'recognised')
        if (beta != None):
            if (self.typeMotion == 'asInput'):
                # Save past control surface angle and corresponding time
                betaTemp = self.beta
                timeTemp = self.time
                
                # Update time member
                self.time = time
                
                # Check for 0 denominator
                if (self.time == timeTemp):
                    raise Exception('Update time equal to previous time')
                else:
                    # Overwrite current angle and rates
                    self.beta = beta
                    self.betaDot = ((self.beta - betaTemp) / 
                                    (self.time - timeTemp))
            else:
                raise Exception('beta specified and typeMotion not ' +
                                'recognised')
        
class VMUnsteadyInput:
    """@brief Contains data for unsteady run of UVLM.
    
    @param WakeLength Length of wake in chordlengths.
    @param DelS non-dim timestep tau = tU/c, s = omega*c/U.
    @param NumChordLengths Number of chord lengths to travel 
    in prescribed simulation.
    @param VelA_G Velocity of reference frame.
    @param OmegaA_G Initial angular vel of reference frame.
    @param OriginA_G Origin of reference frame in G-frame.
    @param PsiA_G Orientation of reference frame in G-frame."""
    
    def __init__(self, VMOPTS, VMINPUT, WakeLength,
                 DelS, NumChordLengths,
                 VelA_G, OmegaA_G,
                 OriginA_G = np.zeros((3),ct.c_double),
                 PsiA_G = np.zeros((3),ct.c_double)):
        
        self.WakeLength = WakeLength
        self.NumChordLengths = NumChordLengths
        self.VelA_G = VelA_G
        self.OmegaA_G = OmegaA_G
        self.OriginA_G = OriginA_G
        self.PsiA_G = PsiA_G
        
        # Physical time step information.
        try:
            self.NumSteps = NumChordLengths/DelS
        except ZeroDivisionError:
            self.NumSteps = 1
        
        self.DelTime = DelS
        self.FinalTime = NumChordLengths*VMINPUT.c / \
                    np.linalg.norm(VMINPUT.U_infty-VelA_G)
        
        # Check Mstar is high enough to discretize wake.
        if self.DelTime != self.WakeLength*VMINPUT.c/VMOPTS.Mstar.value:
            print("DelTime requested is ",self.DelTime,\
                  "\nDelWakePanel is ",self.WakeLength*VMINPUT.c/VMOPTS.Mstar.value,\
                  "\nChanging Mstar to ", self.WakeLength*VMINPUT.c/self.DelTime)
            VMOPTS.Mstar.value = int(self.WakeLength*VMINPUT.c/self.DelTime)
            
        # Set DelTime for VMOPTS (dimensional)
        VMOPTS.DelTime = ct.c_double(self.DelTime*VMINPUT.c/np.linalg.norm(VMINPUT.U_infty-VelA_G))
        
class Gust:
    """@brief Gust.
    
    @details Gusts are specified in the G-frame (inertial).
    """
    
    def __init__(self, uMag, l, r):
        """@brief Initialise a gust.
        
        @param uMag Reference velocity for the gust.
        @param l Gust length.
        @param r Location of gust on G-frame (inertial) y-axis.
        """
        
        self.uMag = uMag
        self.l = l
        self.r = r
    
    def _vel(self, xGust):
        """@brief Gust velocity.
        
        @param xGust Penetration length into the gust.
        @return vel Velocity at xGust.
        """
        
        vel = np.zeros((3),ct.c_double,'C')
        if xGust < 0.0 or xGust > self.l:
            return vel
        else:
            # 1 - cos gust
            vel[2] = 0.5*self.uMag*(1.0-np.cos((2.0*np.pi*xGust)/(self.l)))
            return vel
    
    def _getxGust(self, zeta):
        """@brief Penetration length at grid points.
        
        @param zeta Aerodynamic grid point.
        @return xGust penetration length.
        """
        
        return (zeta[1] - self.r)
    
    def Vels(self, zeta):
        """@brief Calculate gust velocities at grid points.
        
        @param zeta array of grid coordinates (M+1,N+1,3).
        @return uExt array of collocation point gust velocities (M,N,3).
        """
        
        uExt = np.zeros((zeta.shape[0],zeta.shape[1],3),ct.c_double,'C')
        for i in range(zeta.shape[0]):
            for j in range(zeta.shape[1]):
                uExt[i,j,:] = self._vel(self._getxGust(zeta[i,j,:]))
            # END for j
        # END for i
        return uExt
        

if __name__ == '__main__':
    c = 1
    b = 8
    U_mag = 25.0
    alpha = 2*np.pi/180.0
    VMINPUT = VMinput(c, b, U_mag, alpha, 0.0, 0.0, 15.0)
    print(VMINPUT.U_infty)
    
    # Testing gust stuff
    GUST = Gust(1.0,1.0,0.0)
    zeta = np.array([[[0.0,0.0,0.0]]])
    print(GUST.Vels(zeta))
    