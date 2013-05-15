'''
Created on 20 Feb 2013

@author: rjs10
'''
import unittest
import SharPySettings as Settings
import numpy as np
import ctypes as ct
from DerivedTypesAero import VMopts, VMinput
from PyAero.UVLM.Solver.UVLM import Run_Cpp_Solver_UVLM, VMUnsteadyInput
from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps
import matplotlib.pyplot as plt

TestDir = Settings.SharPyProjectDir + 'SharPy/src/PyAero/UVLM/' \
           + 'Solver/test/UVLM/'
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
        AELOPTS = AeroelasticOps(ElasticAxis = -0.5,\
                                gForce =  0.0,\
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
        
        # Define 'aeroelastic' parameters"
        AELOPTS = AeroelasticOps(ElasticAxis = -0.5,\
                                gForce =  0.0,\
                                AirDensity = 1.0)
        
        # Run C++ solver."
        CoeffHistory = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS)
        
        # Check forces from Katz and Plotkin method."
        self.assertAlmostEqual(CoeffHistory[-1,3],2*np.pi*VMINPUT.theta,delta=1e-3)
        self.assertAlmostEqual(-CoeffHistory[-1,2],0,delta=1e-5)

class Test_v_AR8_inKatzAndPlotkin(unittest.TestCase):
    """@brief Test unsteady response to that in Katz and Plotkin."""

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
        VMINPUT = VMinput(c, b, U_mag, alpha, theta, 0.0 ,\
                          NumChordLengths)
        
        # Define unsteady solver parameters.
        WakeLength = 10.0
        DeltaS = 1.0/16.0
        VelA_A = np.array([0, 1.0, 0])
        OmegaA_A = np.array([0, 0, 0])
        VMUNST = VMUnsteadyInput(VMOPTS,VMINPUT,\
                                 WakeLength,\
                                 DeltaS,\
                                 NumChordLengths,\
                                 VelA_A, OmegaA_A)
        
        # Define 'aeroelastic' options.
        AELOPTS = AeroelasticOps(ElasticAxis = -0.5,\
                                gForce =  0.0,\
                                AirDensity = 1.0)
        
        # Run C++ solver.
        CoeffHistory = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS)
        
        
        # Check force for Katz and Plotkin method against textbook at final 
        # time.
        CD = np.loadtxt('UVLM/Katz1p492_AR8_CD.dat',ndmin = 2)
        CL = np.loadtxt('UVLM/Katz1p492_AR8_CL.dat',ndmin = 2)
        
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
    """@brief Test response to control surface motions with the results of 
    Theodorsen and Garrick."""

    def setUp(self):
        # Set SharPy output directory and file root.
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'CtrlSurf_' 


    def tearDown(self):
        pass
        
    
    def test_KatzPlotkin(self):
        """Test Katz and Plotkin method unsteady lift and drag."""
        pass
        

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
    #TestRunner.run(suite2) #run a single suite
    TestRunner.run(alltests) #run all tests