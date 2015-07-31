'''
Created on 20 Feb 2013

@author: rjs10
'''
import unittest
import SharPySettings as Settings
import numpy as np
from DerivedTypesAero import VMopts, VMinput, ControlSurf
from PyAero.UVLM.Solver.UVLM import Run_Cpp_Solver_UVLM, VMUnsteadyInput
from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps
import matplotlib.pyplot as plt
from PyAero.UVLM.Utils.Analytical import TheoGarrAerofoil

TestDir = (Settings.SharPyProjectDir + 'output/tests/PyAero/UVLM/' 
           + 'UVLM/')
#/home/rjs10/Documents/MATLAB/UVLM++/Comparisons

class Test_v_TAT(unittest.TestCase):
    """@brief Test final unsteady lift/drag to TAT for high AR problem."""

    def setUp(self):
        # Set SharPy output directory and file root.
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'TAT_' 


    def tearDown(self):
        pass


    def test_KJ(self):
        """Test final unsteady lift/drag to TAT for high AR problem using KJ
        method."""
        
        M = 1
        N = 1
        Mstar = 200
        VMOPTS = VMopts(M,N,True,Mstar,False,True)
        
        # Define wing parameters.
        c = 1
        b = 10000 #semi-span
        
        # Define free-stream conditions.
        U_mag = 0.0
        alpha = 0.0*np.pi/180.0
        theta = 1.0*np.pi/180.0
        VMINPUT = VMinput(c, b, U_mag, alpha, theta)
        
        # Define unsteady simulation parameters."
        WakeLength = 200.0
        DeltaS = 1.0
        NumChordLengths = 200.0
        VelA_A = np.array([0, 1.0, 0])
        OmegaA_A = np.array([0, 0, 0])
        VMUNST = VMUnsteadyInput(VMOPTS,VMINPUT,\
                                 WakeLength,\
                                 DeltaS,\
                                 NumChordLengths,\
                                 VelA_A, OmegaA_A)
        
        # Define 'aeroelastic' options."
        AELOPTS = AeroelasticOps(ElasticAxis = -0.5,InertialAxis=0.0,
                                AirDensity = 1.0)
        
                
        # Run Cpp-solver."
        CoeffHistory = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS)
        
        # Check force from KJ method"
        self.assertAlmostEqual(CoeffHistory[-1,3],
                               2*np.pi*VMINPUT.theta, delta=1e-3)
        self.assertAlmostEqual(-CoeffHistory[-1,2], 0, delta=1e-5)
        
    
    def test_KatzPlotkin(self):
        """Test final unsteady lift/drag to TAT for high AR problem using
        generalised Katz and Plotkin method."""
        
        M = 1
        N = 1
        Mstar = 200
        VMOPTS = VMopts(M,N,True,Mstar,False,False)
        
        # Define wing geometry parameters."
        c = 1
        b = 10000 #semi-span
        
        # Define free-stream conditions."
        U_mag = 0.0
        alpha = 0.0*np.pi/180.0
        theta = 1.0*np.pi/180.0
        VMINPUT = VMinput(c, b, U_mag, alpha, theta)
        
        # Define unsteady simulation parameters.
        WakeLength = 200.0
        DeltaS = 1.0
        NumChordLengths = 200.0
        VelA_A = np.array([0, 1.0, 0])
        OmegaA_A = np.array([0, 0, 0])
        VMUNST = VMUnsteadyInput(VMOPTS,VMINPUT,\
                                 WakeLength,\
                                 DeltaS,\
                                 NumChordLengths,\
                                 VelA_A, OmegaA_A)
        
        # Define 'aeroelastic' parameters
        AELOPTS = AeroelasticOps(ElasticAxis = -0.5,InertialAxis=0.0,
                                AirDensity = 1.0)
        
        # Run C++ solver."
        CoeffHistory = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS)
        
        # Check forces from Katz and Plotkin method."
        self.assertAlmostEqual(CoeffHistory[-1,3],2*np.pi*VMINPUT.theta,delta=1e-3)
        self.assertAlmostEqual(-CoeffHistory[-1,2],0,delta=1e-5)

class Test_v_AR8_inKatzAndPlotkin(unittest.TestCase):
    """Test unsteady response to that in Katz and Plotkin."""

    def setUp(self):
        # Set SharPy output directory and file root.
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'AR8_' 


    def tearDown(self):
        pass
        
    
    def test_KatzPlotkin(self):
        """Test final unsteady lift/drag for AR8 problem."""
        
        M = 4
        N = 13
        Mstar = 160
        VMOPTS = VMopts(M,N,True,Mstar,False,False,True,0.0,False,4)
        
        # Define wing geometry parameters.
        c = 1
        b = 4 #semi-span
        
        # Define free-stream conditions.
        U_mag = 0.0
        alpha = 0.0*np.pi/180.0
        theta = 5.0*np.pi/180.0
        NumChordLengths = 10.0
        VMINPUT = VMinput(c, b, U_mag, alpha, theta, NumChordLengths)
        
        # Define unsteady solver parameters.
        WakeLength = 10.0
        DeltaS = 1.0/16.0
        VelA_A = np.array([0, 1.0, 0])
        OmegaA_A = np.array([0, 0, 0])
        VMUNST = VMUnsteadyInput(VMOPTS,VMINPUT,
                                 WakeLength,
                                 DeltaS,
                                 NumChordLengths,
                                 VelA_A, OmegaA_A)
        
        # Define 'aeroelastic' options.
        AELOPTS = AeroelasticOps(ElasticAxis = -0.5,InertialAxis=0.0,
                                AirDensity = 1.0)
        
        # Run C++ solver.
        CoeffHistory = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS)
        
        
        # Check force for Katz and Plotkin method against textbook at final 
        # time.
        CD = np.loadtxt('testData/Katz1p492_AR8_CD.dat',ndmin = 2)
        CL = np.loadtxt('testData/Katz1p492_AR8_CL.dat',ndmin = 2)
        
        self.assertAlmostEqual(CoeffHistory[-1,3],CL[-1,1],delta=1e-3)
        self.assertAlmostEqual(-CoeffHistory[-1,2],CD[-1,1],delta=1e-3)

        Plotting = False
        if Plotting == True:
            plt.figure(1)
            a = plt.plot(CD[:,0],CD[:,1],'ks')
            b = plt.plot(CoeffHistory[:,0],-CoeffHistory[:,2],'k-')
            plt.xlim((0, 10.0))
            plt.ylim((0, 0.03))
            plt.xlabel('U*t/c [-]')
            plt.ylabel('C_D')
            plt.legend(('Katz and Plotkin (2001)',
                       'SHARPy'), 'upper right', shadow=True,
                       prop={'size':12})
            
            plt.figure(2)
            c = plt.plot(CL[:,0],CL[:,1],'ks')
            d = plt.plot(CoeffHistory[:,0],CoeffHistory[:,3],'k-')
            plt.xlim((0, 10.0))
            plt.ylim((0, 0.6))
            plt.xlabel('U*t/c [-]')
            plt.ylabel('C_L')
            plt.legend(('Katz and Plotkin (2001)',
                       'SHARPy'), 'lower right', shadow=True,
                       prop={'size':12})
            
            lines = a,b,c,d
            plt.setp(lines,linewidth=2,markersize=6)
            
            plt.show()
        # END if Plotting
        
class Test_v_TheoGarrick(unittest.TestCase):
    """@brief Test response to pitch/plunge/control surface oscillations 
    with the results of Theodorsen and Garrick."""

    def setUp(self):
        # Set SharPy output directory and file root.
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'TheoGarr_' 


    def tearDown(self):
        pass
        
    
    def test_ctrlSurf(self):
        """Test K&P then Joukowski methods unsteady lift and drag."""
        
        # Declare UVLM solver options.
        VMOPTS = VMopts(M = 25,
                        N = 2,
                        ImageMethod = True,
                        Mstar = 250,
                        Steady = False,
                        KJMeth = False,
                        NewAIC = True,
                        DelTime = 0.0,
                        Rollup = False,
                        NumCores = 4)
        
        # Define wing geometry parameters.
        c = 1.0
        b = 1000.0
        theta = 0.0*np.pi/180.0
        
        # Define free-stream conditions.
        U_mag = 1.0
        alpha = 0.0*np.pi/180.0
    
        # Define Control Surface
        k = 1.0
        omega = 2*U_mag*k/c
        betaBar = 0.1*np.pi/180.0
        ctrlSurf = ControlSurf(20, 25, 0, 2, 'sin', betaBar, omega)
        
        # Length of simulation.
        nT = 2.0
        NumChordLengths = nT * np.pi / k
        
        VMINPUT = VMinput(c,
                          b,
                          U_mag,
                          alpha,
                          theta,
                          WakeLength = NumChordLengths,
                          ctrlSurf = ctrlSurf)
        
        # Define unsteady solver parameters.
        VMUNST = VMUnsteadyInput(VMOPTS,
                                 VMINPUT,
                                 WakeLength = 10.0,
                                 DelS = 1.0/25.0,
                                 NumChordLengths = NumChordLengths,
                                 VelA_G = np.array([0, 0, 0]),
                                 OmegaA_G = np.array([0, 0, 0]))
        
        # Define 'aeroelastic' options.
        ElasticAxis = -0.5
        AELOPTS = AeroelasticOps(ElasticAxis = ElasticAxis) # Quarter chord.
        
        # Run C++ solver.
        CoeffHistoryKatz = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS)
        
        # Run from Joukowski method.
        VMOPTS.KJMeth.value = True
        CoeffHistoryJouk = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS)
        
        # Set parameters for analytical result.
        alphaBar = 0.0
        hBar = 0.0
        phi0 = 0.0
        phi1 = 0.0
        phi2 = 0.0
        a = ElasticAxis
        e = 0.6
        k = 1.0
        
        Cl, Cd, alpha, beta, hDot, s = TheoGarrAerofoil(alphaBar, betaBar, hBar,
                                                        phi0, phi1, phi2,
                                                        a, e, k, nT+0.1)
        
        # Delete unused variables
        del hDot
        
        # Interpolate analytical result onto simulated time steps.
        ClTheo = np.interp(CoeffHistoryKatz[:,0]*U_mag/(c/2.0), s, Cl)
        CdGarr = np.interp(CoeffHistoryKatz[:,0]*U_mag/(c/2.0), s, Cd)
        
        # Calculate error as function of time in last period.
        kTimeMid = len(CoeffHistoryKatz[:,0])/2
        ClErr = CoeffHistoryKatz[kTimeMid:,3] - ClTheo[kTimeMid:]
        CdErr = -CoeffHistoryKatz[kTimeMid:,2] - CdGarr[kTimeMid:]
        ClErrJouk = CoeffHistoryJouk[kTimeMid:,3] - ClTheo[kTimeMid:]
        CdErrJouk = -CoeffHistoryJouk[kTimeMid:,2] - CdGarr[kTimeMid:]
        
        # RMS error
        rmsCl = np.sqrt(np.mean(pow(ClErr,2.0)))/max(CoeffHistoryKatz[:,3])
        rmsCd = np.sqrt(np.mean(pow(CdErr,2.0)))/max(-CoeffHistoryKatz[:,2])
        rmsClJouk = ( np.sqrt(np.mean(pow(ClErrJouk,2.0)))
                      / max(CoeffHistoryJouk[:,3]) )
        rmsCdJouk = ( np.sqrt(np.mean(pow(CdErrJouk,2.0)))
                      / max(-CoeffHistoryJouk[:,2]) )
        
        # Check RMS is same as on 23/05/13
        self.assertAlmostEqual(rmsCl, 0.0241589106969,6,'RMS error in lift.')
        self.assertAlmostEqual(rmsCd, 0.136290106597,5,'RMS error in drag.')
    
        Plotting = False
        if Plotting == True:
            # Define beta here for plotting purposes.
            betaSim = betaBar * np.sin(omega * CoeffHistoryKatz[:,0])
            
            plt.subplot(2,1,1)
            plt.grid()
            plt.plot(betaSim,-CoeffHistoryKatz[:,2],'ko-',beta,Cd,'k-.')
            plt.xlabel(r'$\beta$')
            plt.ylabel(r'$C_d$')
            plt.subplot(2,1,2)
            plt.grid()
            plt.plot(betaSim,CoeffHistoryKatz[:,3],'ko-',beta,Cl,'k-.')
            plt.xlabel(r'$\beta$')
            plt.ylabel(r'$C_l$')
            
            plt.show()
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    suite1 = (unittest.TestLoader()
              .loadTestsFromTestCase(Test_v_TAT))
    suite2 = (unittest.TestLoader()
              .loadTestsFromTestCase(Test_v_AR8_inKatzAndPlotkin))
    suite3 = (unittest.TestLoader()
              .loadTestsFromTestCase(Test_v_TheoGarrick))
    alltests = unittest.TestSuite([suite1, suite2, suite3])
    TestRunner = unittest.TextTestRunner(verbosity=2) # creates a test runner
    #TestRunner.run(suite3) #run a single suite
    TestRunner.run(alltests) #run all tests