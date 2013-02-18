'''
Created on 13 Feb 2013

@author: rjs10
'''
import unittest
import SharPySettings as Settings
import DerivedTypesAero
import numpy as np
import ctypes as ct
from VLM import Run_Cpp_Solver_VLM

TestDir = Settings.SharPyProjectDir + 'SharPy/src/PyAero/UVLM/' \
           + 'Solver/test/SteadyVLM/'


class Test_SteadyVLM_v_TAT(unittest.TestCase):
    """@brief Test lift and drag for very large aspect ratio problem."""

    def setUp(self):
        "Set SharPy output directory and file root"
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = ''    
        

    def tearDown(self):
        pass


    def test_2D_Lift_and_Drag(self):
        """check lift equals 2*pi*alpha and drag is zero"""
        
        "set up baseline options"
        VMOPTS = DerivedTypesAero.VMopts(10, 1, True, 1, True, False)
        VMINPUT = DerivedTypesAero.VMinput(1.0, 1000.0, 1.0,\
                                           0.1*np.pi/180.0,\
                                           0.0*np.pi/180.0)
        
        "solve"
        Coeffs = Run_Cpp_Solver_VLM(VMOPTS,VMINPUT)
        print(Coeffs)
        
        
        "check lift then drag"
        """Note: increasing alpha causes geometrical nonlinear effect of wake 
        position to cause discrepancy with TAT.
        TODO: discrepancy between PyUVLM and UVLM++ here in 3rd digit lift."""
        self.assertAlmostEqual(Coeffs[2],2*np.pi*VMINPUT.alpha,4)
        self.assertAlmostEqual(-Coeffs[1],0,5)
        

        
class Test_SteadyVLM_v_UVLMpp(unittest.TestCase):
    """@brief Test lift and drag for very large aspect ratio problem."""

    def setUp(self):
        "Set SharPy output directory and file root"
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = ''    
        

    def tearDown(self):
        pass
        
    def test_3D_Lift_and_Drag_Katz_then_KJ(self):
        """check converged lift and drag matches trusted result from UVLM++"""
        
        "set up baseline options"
        VMOPTS = DerivedTypesAero.VMopts(5, 80, True, 1, True, False)
        VMINPUT = DerivedTypesAero.VMinput(1.0, 16.0, 1.0,\
                                           10.0*np.pi/180.0,\
                                           0.0*np.pi/180.0)
        
        "solve"
        Coeffs = Run_Cpp_Solver_VLM(VMOPTS,VMINPUT)
        print(Coeffs)
        
        "check lift then drag"
        self.assertAlmostEqual(Coeffs[2],0.965402,2)
        self.assertAlmostEqual(-Coeffs[1],0.0118845,3)
        
        
        "check KJ method"
        VMOPTS.KJMeth = ct.c_bool(True)
        Coeffs = Run_Cpp_Solver_VLM(VMOPTS,VMINPUT)
        print(Coeffs)
        
        
        "check lift then drag"
        self.assertAlmostEqual(Coeffs[2],0.993501,5)
        self.assertAlmostEqual(-Coeffs[1],0.0113952,5)
        

class Test_SteadyVLM_v_TAT_KJMeth(unittest.TestCase):
    """@brief Test lift and drag for very large aspect ratio problem."""

    def setUp(self):
        "Set SharPy output directory and file root"
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = ''    
        

    def tearDown(self):
        pass
        
    def test_2D_Lift_and_Drag(self):
        """check converged lift and drag."""
        
        "set up baseline options"
        VMOPTS = DerivedTypesAero.VMopts(10, 1, True, 1, True, True)
        VMINPUT = DerivedTypesAero.VMinput(1.0, 1000.0, 1.0,\
                                           0.1*np.pi/180.0,\
                                           0.0*np.pi/180.0)
        
        "solve"
        Coeffs = Run_Cpp_Solver_VLM(VMOPTS,VMINPUT)
        print(Coeffs)
        
        "check lift then drag"
        self.assertAlmostEqual(Coeffs[2],2*np.pi*VMINPUT.alpha,4)
        self.assertAlmostEqual(-Coeffs[1],0.0,5)
        
    
class Test_SteadyVLM_ZetaDot(unittest.TestCase):
    """@brief Test lift and drag for very large aspect ratio problem."""

    def setUp(self):
        "Set SharPy output directory and file root"
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = ''    
        

    def tearDown(self):
        pass
        
    def test_2D_Lift_and_Drag(self):
        """check lift and drag."""
        
        "set up baseline options"
        VMOPTS = DerivedTypesAero.VMopts(10, 1, True, 1, True, False)
        VMINPUT = DerivedTypesAero.VMinput(1.0, 1000.0, 0.5,\
                                           0.0*np.pi/180.0,\
                                           0.1*np.pi/180.0,\
                                           ZetaDotTest = 0.5)
        
        "solve"
        Coeffs = Run_Cpp_Solver_VLM(VMOPTS,VMINPUT)
        print(Coeffs)
        
        "check lift then drag"
        self.assertAlmostEqual(Coeffs[2],2*np.pi*VMINPUT.theta,4)
        self.assertAlmostEqual(-Coeffs[1],0.0,5)
        
        
        "set up with KJ meth"
        VMOPTS.KJMeth = ct.c_bool(True)
        
        "solve"
        Coeffs = Run_Cpp_Solver_VLM(VMOPTS,VMINPUT)
        print(Coeffs)
        
        "check lift then drag"
        self.assertAlmostEqual(Coeffs[2],2*np.pi*VMINPUT.theta,4)
        self.assertAlmostEqual(-Coeffs[1],0.0,5)
        
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    suite1 = unittest.TestLoader().loadTestsFromTestCase(\
                Test_SteadyVLM_v_TAT)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(\
                Test_SteadyVLM_v_UVLMpp)
    suite3 = unittest.TestLoader().loadTestsFromTestCase(\
                Test_SteadyVLM_v_TAT_KJMeth)
    suite4 = unittest.TestLoader().loadTestsFromTestCase(\
                Test_SteadyVLM_ZetaDot)
    alltests = unittest.TestSuite([suite1, suite2, suite3, suite4])
    #alltests.run(unittest.defaultTestResult())
    TestRunner = unittest.TextTestRunner(verbosity=2) # creates a test runner
    #TestRunner.run(suite3) #run a single suite
    TestRunner.run(alltests) #run all tests