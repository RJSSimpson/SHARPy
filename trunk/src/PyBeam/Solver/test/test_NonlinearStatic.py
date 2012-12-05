'''@package PyBeam.Solver.test.test_NonlinearStatic
@brief      Python unittest test cases for nonlinear static solver.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       23/11/2012
@pre        None
@warning    None'''
import unittest
import os
import SharPySettings as Settings
import DerivedTypes

TestDir = Settings.SharPyProjectDir + 'SharPy/src/PyBeam/' \
           + 'Solver/test/NonlinearStatic/'

class TestNonlinearStatic_v_Executable(unittest.TestCase):
    """@brief Test Python NonlinearStatic_F90 v F90 Executable for 'TPY0' case
    defined in input_rob.f90."""

    def setUp(self):
        """@brief Clean and compile F90 executable, and setup PyBeam (including
        clean and compile shared library), to run/emulate input_rob.f90 'TPY0' 
        case.
        
        @pre the test case in input_rob.f90 must be set to 'TPY0' manually.
        TODO: Somehow automate the compilation of F90 exec with 'TYP0' case."""
        
        "Clean and compile .so shared library using existing Makefiles"
        os.chdir(Settings.SharPyProjectDir + 'BeamLib/src/wrapper/')
        os.system('make clean')
        os.system('make all')
        
        "Clean and compile F90 executable using existing Makefiles"
        os.chdir(Settings.SharPyProjectDir + \
                 'BeamLib/src/fortran/install/linux')
        os.system('make clean')
        os.system('make all')
        
        "cd to F90 test results directory and run executable"
        os.chdir(TestDir)
        os.system('rm *_SOL*') #remove existing solution files, if any
        ExecDir = Settings.SharPyProjectDir + 'BeamLib/src/fortran/bin/'
        os.system('%s./xbeam.exe' % (ExecDir))
        
        "Set SharPy output directory and file root"
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'PyBeam'


    def tearDown(self):
        """@brief Clean F90 executable compiled for 'TPY0' case."""
        os.chdir(Settings.SharPyProjectDir + \
                 'BeamLib/src/fortran/install/linux')
        os.system('make clean')


    def test_TipDispRot(self):
        """@brief Tests PyBeam tip displacement against F90 Exec tip 
        displacement."""
        
        import NonlinearStatic # imported after clean/make process
        
        "Test if correctly named solution file is created"
        self.assertTrue(os.path.isfile('TPY0_SOL112_def.txt'),'TPY0 solution'\
                        +' does not exist: check testcase in input_rob.f90'\
                        +' is set to \'TPY0\' ')
        
        "Read executable output"
        f = open('TPY0_SOL112_def.txt')
        fLines = f.readlines() #save all lines
        f.close()
        fLast = str(fLines[-1]) #isolate final line as string
        fWords = fLast.split() #split into words
        
        """Set up Xbopts for nonlinear static analysis defined in input_rob.f90
        TPY0 test case"""
        XBOPTS = DerivedTypes.Xbopts()
        XBOPTS.Solution.value = 112 
        XBOPTS.NumLoadSteps.value = 10
        XBOPTS.MinDelta.value = 1e-04       
        """Set up Xbinput for nonlinear static analysis defined in input_rob.f90
        TPY0 test case"""
        XBINPUT = DerivedTypes.Xbinput()
        XBINPUT.NumElems = 8
        XBINPUT.BeamLength = 16.0
        XBINPUT.BeamStiffness[0,0] = 1.0e+09
        XBINPUT.BeamStiffness[1,1] = 1.0e+09
        XBINPUT.BeamStiffness[2,2] = 1.0e+09
        XBINPUT.BeamStiffness[3,3] = 1.0e+04
        XBINPUT.BeamStiffness[4,4] = 2.0e+04
        XBINPUT.BeamStiffness[5,5] = 4.0e+06
        XBINPUT.BeamMass[0,0] = 0.75
        XBINPUT.BeamMass[1,1] = 0.75
        XBINPUT.BeamMass[2,2] = 0.75
        XBINPUT.BeamMass[3,3] = 0.1
        XBINPUT.BeamMass[4,4] = 0.001
        XBINPUT.BeamMass[5,5] = 0.001
        XBINPUT.ForceStatic[2] = 800
        
        "Solve using F90 subroutine"
        NonlinearStatic.Solve_F90(XBINPUT,XBOPTS)
        
        "Read PyBeam output"
        g = open('PyBeam_SOL112_def.dat')
        gLines = g.readlines() #save all lines
        g.close()
        gLast = str(gLines[-1]) #isolate final line as string
        gWords = gLast.split() #split into words
        
        "Compare in-plane displacements and rotations"
        self.assertTrue(fWords[2].lower() ==gWords[2], 'In-plane x-'\
                        + 'displacement does not match')
        self.assertTrue(fWords[4].lower() ==gWords[4], 'In-plane z-'\
                        + 'displacement does not match')
        self.assertTrue(fWords[6].lower() ==gWords[6], 'In-plane y-'\
                        + 'rotation does not match')


class TestNonlinearStatic_v_GeradinCardonna(unittest.TestCase):
    """@brief Test Python NonlinearStatic_F90 against Geradin and Cardonna's
    test case used by Dr. Rafa Palacios."""
    
    def setUp(self):
        """@brief cd to and set SharPy output directory to TestDir,
         and set file root to PyBeam2"""
        os.chdir(TestDir)  
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'PyBeam2'
    

    def tearDown(self):
        pass
    
    
    def test_2and3nodedElems(self):
        """@brief Test results of PyBeam for 2- and 3- noded elements"""
        
        import NonlinearStatic # imported after clean/make process
        
        """Set up Xbopts for nonlinear static analysis defined in 
        NonlinearStatic/testcases.pdf case 1.1."""
        XBOPTS = DerivedTypes.Xbopts()
        XBOPTS.Solution.value = 112 
        XBOPTS.NumLoadSteps.value = 10
        XBOPTS.MinDelta.value = 1e-05       
        """Set up Xbinput for nonlinear static analysis defined in 
        NonlinearStatic/testcases.pdf case 1.1."""
        XBINPUT = DerivedTypes.Xbinput()
        XBINPUT.NumElems = 20
        XBINPUT.BeamLength = 5.0
        XBINPUT.BeamStiffness[0,0] = 4.8e+08/XBINPUT.BeamLength
        XBINPUT.BeamStiffness[1,1] = 3.231e+08/XBINPUT.BeamLength
        XBINPUT.BeamStiffness[2,2] = 3.231e+08/XBINPUT.BeamLength
        XBINPUT.BeamStiffness[3,3] = 1.0e+06/XBINPUT.BeamLength
        XBINPUT.BeamStiffness[4,4] = 9.346e+06/XBINPUT.BeamLength
        XBINPUT.BeamStiffness[5,5] = 9.346e+06/XBINPUT.BeamLength
        XBINPUT.BeamMass[0,0] = 100
        XBINPUT.BeamMass[1,1] = 100
        XBINPUT.BeamMass[2,2] = 100
        XBINPUT.BeamMass[3,3] = 10
        XBINPUT.BeamMass[4,4] = 0.0 #Neglect the cross-section bending inertia
        XBINPUT.BeamMass[5,5] = 0.0 #Neglect the cross-section bending inertia
        XBINPUT.ForceStatic[2] = 6e+05
        
        NonlinearStatic.Solve_F90(XBINPUT, XBOPTS)
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    #unittest.main(verbosity=2)
    suite1 = unittest.TestLoader().loadTestsFromTestCase(\
                TestNonlinearStatic_v_Executable)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(\
                TestNonlinearStatic_v_GeradinCardonna)
    alltests = unittest.TestSuite([suite1, suite2])
    #alltests.run(unittest.defaultTestResult())
    TestRunner = unittest.TextTestRunner(verbosity=2) # creates a test runner
    TestRunner.run(suite1) #run a single suite
    #TestRunner.run(alltests) #run all tests