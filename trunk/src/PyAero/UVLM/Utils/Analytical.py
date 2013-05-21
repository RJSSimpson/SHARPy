'''PyAero.UVLM.Utils.Analytical
@brief      Analytical solution for comparison with UVLM.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       21/05/2013
@pre        None
@warning    None
'''
from scipy.constants.constants import pi
import numpy as np
from numpy import sin, cos, arccos, arctan2
from scipy.special import kv
import matplotlib.pyplot as plt
from numpy.lib.type_check import real, imag

def C(k):
    """@brief Theodorsen's function.
    @param k Reduced frequency.
    @return C Function value, complex type.
    """
    
    # Theodorsen's Function
    if k < 0.0:
        raise Exception('Reduced frequency should not be negative.')
    elif k == 0.0:
        return complex(1.0,0.0)
    else:
        return kv(1, complex(0, k)) / (kv(0, complex(0, k)) + 
                                       kv(1, complex(0, k)))
    

def TheoGarrAerofoil(alphaBar,betaBar,hBar,phi0,phi1,phi2,a,e,k,nT):
    """@brief Lift and drag of thin aerofoil oscillating in pitch, plunge and
    flap angle.
    @param alphaBar amplitude of pitching oscillations.
    @param betaBar amplitude of flap oscillations.
    @param hBar amplitude of plunging oscillations.
    @param phi0 Phase angle in pitch.
    @param phi1 Phase angle in flap angle.
    @param phi2 Phase angle in plunge.
    @param a Non-dim offset of pitching axis aft of mid-chord.
    @param e Non-dim offset of flap hinge aft of mid-chord.
    @param k Reduced frequency, k = omega*b/U.
    @param nT Number of periods in output."""
    
    # Define dummy variables for dimensionality.
    U = 1.0
    b = 1.0
    
    # Time parameters.
    omega = k * U / b
    tFinal = 2 * pi * nT / omega
    t = np.linspace(0.0, tFinal, nT * 100, True)
    
    # Kinematic Parameters
    # Sinusoidal variations, take imaginary part of Theodorsen's equations.
    alpha = alphaBar * sin(omega*t + phi0)
    alphaDot = omega * alphaBar * cos(omega*t + phi0)
    alphaDotDot = - pow(omega,2.0) * alphaBar * sin(omega*t + phi0)
    h = hBar * sin(omega*t + phi1)
    hDot = omega * hBar * cos(omega*t + phi1)
    hDotDot = - pow(omega,2.0) * hBar * sin(omega*t + phi1)
    beta = betaBar * sin(omega*t + phi2)
    betaDot = omega * betaBar * cos(omega*t + phi2)
    betaDotDot = - pow(omega,2.0) * betaBar * sin(omega*t + phi2)
    
    # Theodorsen's constants
    F1 = e * arccos(e) - 1.0/3.0 * (2 + pow(e, 2.0)) * np.sqrt(1 - pow(e, 2.0))
    F2 = (e * (1 - pow(e, 2.0)) - 
          (1 + pow(e, 2.0)) * np.sqrt(1 - pow(e, 2.0)) * arccos(e) +
          e*pow(arccos(e), 2.0))
    F3 = (- (1.0/8.0 + pow(e, 2.0)) * pow(arccos(e),2.0) +
            1.0/4.0 * e * np.sqrt(1 - pow(e, 2.0)) * arccos(e)
                    * (7 + 2*pow(e, 2.0))
          - 1.0/8.0 * (1 - pow(e, 2.0)) * (5 * pow(e, 2.0) + 4))
    F4 = e * np.sqrt(1 - pow(e, 2.0)) - arccos(e)
    F5 = (-(1 - pow(e, 2.0)) - pow(arccos(e),2.0) 
           + 2 * e * np.sqrt(1 - pow(e, 2.0)) * arccos(e))
    F6 = F2
    F7 = (-(1.0/8.0 + pow(e, 2.0)) * arccos(e) 
          + 1.0/8.0 * e * np.sqrt(1 - pow(e, 2.0)) * (7 + 2 * pow(e, 2.0)))
    F8 = (-1.0/3.0 * np.sqrt(1 - pow(e, 2.0)) * (2 * pow(e, 2.0) + 1) 
          + e * arccos(e))
    F9 = 1.0/2.0 * (1.0/3.0 * pow(np.sqrt(1 - pow(e, 2.0)),3) + a * F4)
    F10 = np.sqrt(1 - pow(e, 2.0)) + arccos(e)
    F11 = (1 - 2 * e) * arccos(e) + (2 - e) * np.sqrt(1 - pow(e, 2.0))
    F12 = np.sqrt(1 - pow(e, 2.0)) * (2 + e) - arccos(e) * (2 * e + 1)
    F13 = 1.0/2.0 * (-F7 - (e - a) * F1)
    # F14 = -arccos(e) + e*np.sqrt(1-pow(e, 2.0)) 
    F14 = 1.0/16.0 + 1.0/2.0 * a * e
    
    F6 = F2

    # Garrick's constants
    F15 = F4 + F10
    F16 = F1 - F8 - (e - a) * F4 + 1.0/2.0 * F11
    F17 = -2 * F9 - F1 + (a - 0.5) * F4
    F18 = F5 - F4 * F10
    F19 = F4 * F11
    F20 = -np.sqrt(1 - pow(e, 2.0)) + arccos(e)
    
    M = (  2 * real(C(k)) * (U * alphaBar * cos(phi0) - hBar * omega * sin(phi2)
                             - b * (0.5 - a) * alphaBar * omega * sin(phi0))
         - 2 * imag(C(k)) * (U * alphaBar * sin(phi0) + hBar * omega * cos(phi2)
                             + b * (0.5 - a) * alphaBar * omega * cos(phi0))
         + b * alphaBar * omega * sin(phi0)  )
    
    N = (  2 * real(C(k)) * (U * alphaBar * sin(phi0) + hBar * omega * cos(phi2)
                             + b * (0.5 - a) * alphaBar * omega * cos(phi0))
         + 2 * imag(C(k)) * (U * alphaBar * cos(phi0) - hBar * omega * sin(phi2)
                             - b * (0.5 - a) * alphaBar * omega * sin(phi0))
         - b * alphaBar * omega * cos(phi0)  )
    
    # Loop through time steps
    p = ( -(pi*alphaDot + pi*hDotDot - pi*a*alphaDotDot 
            - F4*betaDot - F1*betaDotDot)
          - 2*pi*real(C(k)) * (alpha + hDot + (0.5 - a)*alphaDot + F10/pi*beta
                               + F11/(2*pi)*betaDot) 
          - 2*pi*imag(C(k))/omega * (alphaDot + hDotDot + alphaDotDot 
                                     + F10/pi*betaDot + F11/(2*pi)*betaDotDot) )
    
    pBeta = ( -(-F4*alphaDot - F4*hDotDot + F9*alphaDotDot - F5/(2*pi)*betaDot
                - F2/(2*pi)*betaDotDot)
              - 2*np.sqrt(1 - pow(e,2.0)) * (0.5*(1 - e)*alphaDot
                                             + np.sqrt(1 - pow(e,2.0))/pi*beta
                                             + F10/(2*pi)*(1 - e)*betaDot)
              - 2*F20*real(C(k))*(alpha + hDot + (0.5 - a)*alphaDot
                                  + F10/pi*beta + F11/(2*pi)*betaDot)
              - 2*F20*imag(C(k))/omega*(alphaDot + hDotDot
                                        + (0.5 - a)*alphaDotDot
                                        + F10/pi*betaDot
                                        + F11/(2*pi)*betaDotDot) )
    
    S = (np.sqrt(2.0 * b) / 2.0) * (M * sin(omega * t) + 
                                    N * cos(omega * t)  )

    d = - alpha*p - pi*pow(S, 2.0) - beta*pBeta
    
    Cl = - p/(pow(U,2.0)*b)
    Cd = d/(pow(U,2.0)*b)
        
    return Cl, Cd, alpha, beta, hDot


if __name__ == '__main__':
    # Test Theodorsen function.
    if False:
        k = np.linspace(0,1,101,True)
        Ck = np.zeros_like(k,complex)
        for kIter in range(k.shape[0]):
            Ck[kIter] = C(k[kIter])
        plt.plot(k,real(Ck),k,imag(Ck))
        plt.grid()
        plt.xlabel('k')
        plt.legend(['Re(C(k))','Im(C(k))'])
        plt.show()
    
    # Test TheoGarr.
    if True:
        alphaBar = 0.0*pi/180.0
        betaBar = 0.1*pi/180.0
        hBar = 0.0
        phi0 =0.0*pi/180.0
        phi1 = 0.0*pi/180.0
        phi2 = 0.0*pi/180.0
        a = -0.5
        e = 0.6
        k = 0.1
        nT = 1.0
        Cl, Cd, alpha, beta, hDot = TheoGarrAerofoil(alphaBar, betaBar, hBar,
                                                     phi0, phi1, phi2,
                                                     a, e, k, nT)
        # Read old UVLM results.
        Katz = np.loadtxt('/home/rjs10/Documents/MATLAB/Theodorsen/' + 
                          'FlappedAerofoil/KatzCoeff50.dat',ndmin = 2)
        Jouk = np.loadtxt('/home/rjs10/Documents/MATLAB/Theodorsen/' + 
                          'FlappedAerofoil/JoukCoeff50.dat',ndmin = 2)
        
        # Plot simulation against analytical.
        betaSim = betaBar*sin(2*k*Katz[:,0])
        plt.plot(betaSim*180.0/pi,Katz[:,1],betaSim*180.0/pi,Jouk[:,1],
                 beta*180.0/pi,Cd,'k+')
        plt.show()
        
        plt.plot(betaSim*180.0/pi,Katz[:,3],betaSim*180.0/pi,Jouk[:,3],
                 beta*180.0/pi,Cl,'k+')
        plt.show()
        
#        alpha_eff = arctan2(hDot,1.0)
#        plt.plot(beta*180.0/pi, Cl)
#        plt.grid()
#        plt.xlabel(r'$\alpha$')
#        plt.ylabel(r'$C_l$')
#        plt.show()
#        
#        plt.plot(beta*180.0/pi, Cd)
#        plt.grid()
#        plt.xlabel(r'$\alpha$')
#        plt.ylabel(r'$C_d$')
#        plt.show()
        