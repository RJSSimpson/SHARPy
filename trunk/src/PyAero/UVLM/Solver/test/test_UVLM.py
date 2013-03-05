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
        "Set SharPy output directory and file root"
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'TAT_' 


    def tearDown(self):
        pass


    def test_KJ(self):
        """Test final unsteady lift/drag to TAT for high AR problem."""
        M = 1
        N = 1
        Mstar = 200
        VMOPTS = VMopts(M,N,True,Mstar,False,True)
        
        "define wing geom"
        c = 1
        b = 10000 #semi-span
        
        "free stream conditions"
        U_mag = 0.0
        alpha = 0.0*np.pi/180.0
        theta = 1.0*np.pi/180.0
        VMINPUT = VMinput(c, b, U_mag, alpha, theta)
        
        "unsteady stuff"
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
        
        "Aeroelastic options"
        AELOPTS = AeroelasticOps(ElasticAxis = -0.5,\
                                gForce =  0.0,\
                                AirDensity = 1.0)
        
        "run solver"
        CoeffHistory = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS)
        
        "Check KJ method"
        self.assertAlmostEqual(CoeffHistory[-1,3],2*np.pi*VMINPUT.theta,delta=1e-3)
        self.assertAlmostEqual(-CoeffHistory[-1,2],0,delta=1e-5)
        
    
    def test_KatzPlotkin(self):
        """Test final unsteady lift/drag to TAT for high AR problem."""
        M = 1
        N = 1
        Mstar = 200
        VMOPTS = VMopts(M,N,True,Mstar,False,False)
        
        "define wing geom"
        c = 1
        b = 10000 #semi-span
        
        "free stream conditions"
        U_mag = 0.0
        alpha = 0.0*np.pi/180.0
        theta = 1.0*np.pi/180.0
        VMINPUT = VMinput(c, b, U_mag, alpha, theta)
        
        "unsteady stuff"
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
        
        "Aeroelastic options"
        AELOPTS = AeroelasticOps(ElasticAxis = -0.5,\
                                gForce =  0.0,\
                                AirDensity = 1.0)
        
        "run solver"
        CoeffHistory = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS)
        
        "Check Katz and Plotkin method"
        self.assertAlmostEqual(CoeffHistory[-1,3],2*np.pi*VMINPUT.theta,delta=1e-3)
        self.assertAlmostEqual(-CoeffHistory[-1,2],0,delta=1e-5)

class Test_v_AR8_inKatzAndPlotkin(unittest.TestCase):
    """@brief Test unsteady response to that in Katz and Plotkin."""

    def setUp(self):
        "Set SharPy output directory and file root"
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
        
        "define wing geom"
        c = 1
        b = 4 #semi-span
        
        "free stream conditions"
        U_mag = 0.0
        alpha = 0.0*np.pi/180.0
        theta = 5.0*np.pi/180.0
        NumChordLengths = 10.0
        VMINPUT = VMinput(c, b, U_mag, alpha, theta, 0.0 ,\
                          NumChordLengths)
        
        "unsteady stuff"
        WakeLength = 10.0
        DeltaS = 1.0/16.0
        VelA_A = np.array([0, 1.0, 0])
        OmegaA_A = np.array([0, 0, 0])
        VMUNST = VMUnsteadyInput(VMOPTS,VMINPUT,\
                                 WakeLength,\
                                 DeltaS,\
                                 NumChordLengths,\
                                 VelA_A, OmegaA_A)
        
        "Aeroelastic options"
        AELOPTS = AeroelasticOps(ElasticAxis = -0.5,\
                                gForce =  0.0,\
                                AirDensity = 1.0)
        
        "run solver"
        CoeffHistory = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS)
        
        
        "Check Katz and Plotkin method"
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
            plt.legend(('Katz and Plotkin (2001)',\
                        'SHARPy'),\
               'upper right',shadow=True,prop={'size':12})
            
            plt.figure(2)
            c = plt.plot(CL[:,0],CL[:,1],'ks')
            d = plt.plot(CoeffHistory[:,0],CoeffHistory[:,3],'k-')
            plt.xlim((0, 10.0))
            plt.ylim((0, 0.6))
            plt.xlabel('U*t/c [-]')
            plt.ylabel('C_L')
            plt.legend(('Katz and Plotkin (2001)',\
                        'SHARPy'),\
               'lower right',shadow=True,prop={'size':12})
            
            lines = a,b,c,d
            plt.setp(lines,linewidth=2,markersize=6)
            
            plt.show()
        # END if Plotting
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    suite1 = unittest.TestLoader().loadTestsFromTestCase(\
                Test_v_TAT)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(\
                Test_v_AR8_inKatzAndPlotkin)
    alltests = unittest.TestSuite([suite1, suite2])
    TestRunner = unittest.TextTestRunner(verbosity=2) # creates a test runner
    #TestRunner.run(suite2) #run a single suite
    TestRunner.run(alltests) #run all tests