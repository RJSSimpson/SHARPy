'''@package PyBeam.Solver.test.test_NonlinearDynamic
@brief      Python unittest test cases for nonlinear dynamic solver.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       23/11/2012
@pre        None
@warning    None'''
import unittest
import os
import sys
import SharPySettings as Settings
import DerivedTypes

TestDir = Settings.SharPyProjectDir + 'SharPy/src/PyBeam/' \
           + 'Solver/test/NonlinearDynamic/'
           
class TestNonlinearStatic_v_Executable(unittest.TestCase):
        """@brief Test Python NonlinearDynamic_F90 v F90 Executable for 'TPY0' case
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
            sys.stderr.write('test_TipDispRot '\
                             + '(__main__.TestNonlinearDynamic_v_Executable)\n')
            
            
        def test_TipDispRot(self):
            """@brief Tests PyBeam tip displacement against F90 Exec tip 
            displacement."""
            
            import NonlinearDynamic # imported after clean/make process
        
            ExecOutFile = 'TPY0_SOL312_dyn.txt' #same run used for dynamic test
            
            "Test if correctly named solution file is created"
            self.assertTrue(os.path.isfile(ExecOutFile),'TPY0 solution'\
                            +' does not exist: check testcase in input_rob.f90'\
                            +' is set to \'TPY0\' ')
        
            "Read executable output"
            f = open(ExecOutFile)
            fLines = f.readlines() #save all lines
            f.close()
            fLast = str(fLines[-1]) #isolate final line as string
            fWords = fLast.split() #split into words
            
            """Set up Xbopts for nonlinear dynamic analysis defined in input_rob.f90
            TPY0 test case"""
            XBOPTS = DerivedTypes.Xbopts()
            XBOPTS.Solution.value = 312 
            XBOPTS.NumLoadSteps.value = 10
            XBOPTS.MinDelta.value = 1e-04   
            XBOPTS.NumGauss.value = 2
            XBOPTS.NewmarkDamp.value = 0.01
            XBOPTS.PrintInfo.value = False
             
            """Set up Xbinput for nonlinear analysis defined in input_rob.f90
            TPY0 test case"""
            XBINPUT = DerivedTypes.Xbinput(3,8)
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
            XBINPUT.ForceStatic[-1,2] = 80
            
            "dynamic parameters"
            XBINPUT.t0 = 0.0
            XBINPUT.tfin = 1.0
            XBINPUT.dt = 0.01
            XBINPUT.Omega = 20.0
            XBINPUT.ForceDyn[-1,2] = 160
            
            
            "Solve using F90 subroutine"
            NonlinearDynamic.Solve_F90(XBINPUT,XBOPTS)
            
            
            "Read PyBeam output"
            g = open('PyBeam_SOL312_dyn.dat')
            gLines = g.readlines() #save all lines
            g.close()
            gLast = str(gLines[-1]) #isolate final line as string
            gWords = gLast.split() #split into words
            
            
            "Compare in-plane displacements and rotations with "
            self.assertAlmostEqual(float(fWords[0]), float(gWords[0]), 3,\
                                    'Final time does not match')
            self.assertAlmostEqual(float(fWords[1]), float(gWords[1]), 3,\
                                    'Final x-displacement does not match')
            self.assertAlmostEqual(float(fWords[3]), float(gWords[3]), 3,\
                                    'Final y-displacement does not match')
            self.assertAlmostEqual(float(fWords[5]), float(gWords[5]), 3,\
                                    'Final z-rotation does not match')

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()