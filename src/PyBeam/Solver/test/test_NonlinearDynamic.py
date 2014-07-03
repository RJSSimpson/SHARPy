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
import numpy as np
import ctypes as ct
import matplotlib.pyplot as plt

TestDir = Settings.SharPyProjectDir + 'output/tests/PyBeam/NlnDynamic/'
           
class TestNonlinearDynamic_v_Executable(unittest.TestCase):
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
        XBOPTS.PrintInfo.value = True #False
         
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
            
class TestNonlinearDynamic_v_PalaciosCesnik(unittest.TestCase):
    """@brief Test nonlinear dynamic solution against that from Palacios, R. 
    and Cesnik, C. E. S. (2009) Structural Models for Flight Dynamics Analysis
     of VFA. 50th AIAA/ASME/ASCE/AHS/ASC Structures, SD and Mat. Conf. 
     AIAA 2009-2403. Palm Springs, CA, USA. 4 - 7th May."""
     
    def setUp(self):
        "Set SharPy output directory and file root"
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'PyBeam2'
        
    
    def tearDown(self):
        pass
    
    
    def test_DynTipDispRot(self):
        """@brief Compare values of F90 to those in Palacios and Cesnik."""
        
        Settings.OutputFileRoot = 'PyBeam2_F90'
        import NonlinearDynamic # imported after clean/make process
         
        """Set up Xbopts for nonlinear static analysis defined in 
        NonlinearStatic/testcases.pdf case 1.2."""
        XBOPTS = DerivedTypes.Xbopts()
        XBOPTS.Solution.value = 312 
        XBOPTS.NumLoadSteps.value = 10
        XBOPTS.MinDelta.value = 1e-05
        XBOPTS.FollowerForce.value = False
        XBOPTS.NumGauss.value = 1
                
        """Set up Xbinput for nonlinear static analysis defined in 
        NonlinearStatic/testcases.pdf case 1.2."""
        XBINPUT = DerivedTypes.Xbinput(3,10)
        XBINPUT.BeamLength = 5.0
        XBINPUT.BeamStiffness[0,0] = 4.8e+08
        XBINPUT.BeamStiffness[1,1] = 3.231e+08
        XBINPUT.BeamStiffness[2,2] = 3.231e+08
        XBINPUT.BeamStiffness[3,3] = 1.0e+06
        XBINPUT.BeamStiffness[4,4] = 9.346e+06
        XBINPUT.BeamStiffness[5,5] = 9.346e+06
        XBINPUT.BeamMass[0,0] = 100
        XBINPUT.BeamMass[1,1] = 100
        XBINPUT.BeamMass[2,2] = 100
        XBINPUT.BeamMass[3,3] = 10
        XBINPUT.BeamMass[4,4] = 0.0001
        XBINPUT.BeamMass[5,5] = 0.0001
         
        "Dynamic parameters"
        XBINPUT.t0 = 0.0
        XBINPUT.tfin = 2.0
        XBINPUT.dt = 0.001
        XBINPUT.Omega = 20.0
        XBINPUT.ForceDyn[-1,2] = 60e+03
        XBINPUT.ForcingType = 'RampSin'
        XBINPUT.RampTime = 1.0
     
        XBOPTS.PrintInfo.value = False
        NonlinearDynamic.Solve_F90(XBINPUT,XBOPTS)
 
        "Plot non-zero tip displacements v Time from file"
 
        "Read from file"
        Dyn = np.loadtxt(Settings.OutputDir + Settings.OutputFileRoot + '_SOL'+\
                        str(XBOPTS.Solution.value) + '_dyn.dat')
     
         
        "Assert z-displacement at 1, 1.5 & 2 secs"
        """Without the original data set the precision of the Figure data is
        limited by jpeg resolution and human error."""
        FigureData = np.loadtxt('testData/PalaciosCesnikFigure_R3.dat')
        self.assertAlmostEqual(Dyn[1000,3],FigureData[0,1], 1,\
                                'z-displacement does not match')
        self.assertAlmostEqual(Dyn[1500,3],FigureData[1,1], 1,\
                                'z-displacement does not match')
        self.assertAlmostEqual(Dyn[2000,3],FigureData[2,1], 1,\
                                'z-displacement does not match')
         
        PlotThings = True
        if PlotThings == True:
            "plot displacements"
            plt.figure(1)
            plt.subplot(311)
            plt.plot(Dyn[:,0],Dyn[:,1])
            plt.grid()
            plt.ylabel('R_1')
            plt.ylim(4.95,5.05)
            plt.subplot(312)
            plt.plot(Dyn[:,0],Dyn[:,3])
            plt.grid()
            plt.ylabel('R_3')
            plt.subplot(313)
            plt.plot(Dyn[:,0],Dyn[:,5])
            plt.grid()
            plt.ylabel('\phi_2')
            plt.xlabel('time')
            #plt.show()
             
            "Run Steps Calc"
            NonlinearDynamic.Solve_F90_steps(XBINPUT,XBOPTS)
            "Read from file"
            Dyn = np.loadtxt(Settings.OutputDir + Settings.OutputFileRoot + '_SOL'+\
                        str(XBOPTS.Solution.value) + '_dyn.dat')
            "plot displacements"
            #plt.figure(1)
            plt.subplot(311)
            plt.plot(Dyn[:,0],Dyn[:,1],'g--')
            #plt.grid()
            plt.ylabel('R_1')
            plt.ylim(4.95,5.05)
            plt.subplot(312)
            plt.plot(Dyn[:,0],Dyn[:,3],'g--')
            #plt.grid()
            plt.ylabel('R_3')
            plt.subplot(313)
            plt.plot(Dyn[:,0],Dyn[:,5],'g--')
            #plt.grid()
            plt.ylabel('\phi_2')
            plt.xlabel('time')
             
             
            "Run Py Calc"
            NonlinearDynamic.Solve_Py(XBINPUT,XBOPTS)
            "Read from file"
            Dyn = np.loadtxt(Settings.OutputDir + Settings.OutputFileRoot + '_SOL'+\
                        str(XBOPTS.Solution.value) + '_dyn.dat')
            "plot displacements"
            #plt.figure(1)
            plt.subplot(311)
            plt.plot(Dyn[:,0],Dyn[:,1],'r.')
            plt.legend(('.f90 solver','.f90 w/ steps in Python', 'Python'),\
                        'upper center', shadow=True)
            #plt.grid()
            plt.ylabel('R_1')
            plt.ylim(4.95,5.05)
            plt.subplot(312)
            plt.plot(Dyn[:,0],Dyn[:,3],'r.')
            #plt.grid()
            plt.ylabel('R_3')
            plt.subplot(313)
            plt.plot(Dyn[:,0],Dyn[:,5],'r.')
            #plt.grid()
            plt.ylabel('\phi_2')
            plt.xlabel('time')
            plt.show()
             
            "plot force-amplitude-in-time"
            plt.figure(2)
            Force = np.loadtxt(Settings.OutputDir + Settings.OutputFileRoot + \
                               '_SOL'+ str(XBOPTS.Solution.value) +\
                               '_force.dat', skiprows=3)
            plt.plot(Force[:,0],Force[:,1])
            plt.grid()
            plt.ylabel('F_{time}')
            plt.show()
 
             
    def test_DynTipDispRot_F90_steps(self):
        """@brief Compare to values of F90 steps to those of F90."""
          
        "Read file from previous test (F90)"
        DynF90 = np.loadtxt(Settings.OutputDir + Settings.OutputFileRoot+'_F90'+
                            '_SOL' + str(312) + '_dyn.dat')
         
        # Update output filename
        Settings.OutputFileRoot = 'PyBeam2_F90steps'
          
        import NonlinearDynamic
          
        """Set up Xbopts for nonlinear static analysis defined in 
        NonlinearStatic/testcases.pdf case 1.2."""
        XBOPTS = DerivedTypes.Xbopts()
        XBOPTS.Solution.value = 312 
        XBOPTS.NumLoadSteps.value = 10
        XBOPTS.MinDelta.value = 1e-05
        XBOPTS.FollowerForce.value = False
        XBOPTS.NumGauss.value = 1
                 
        """Set up Xbinput for nonlinear static analysis defined in 
        NonlinearStatic/testcases.pdf case 1.2."""
        XBINPUT = DerivedTypes.Xbinput(2,20)
        XBINPUT.BeamLength = 5.0
        XBINPUT.BeamStiffness[0,0] = 4.8e+08
        XBINPUT.BeamStiffness[1,1] = 3.231e+08
        XBINPUT.BeamStiffness[2,2] = 3.231e+08
        XBINPUT.BeamStiffness[3,3] = 1.0e+06
        XBINPUT.BeamStiffness[4,4] = 9.346e+06
        XBINPUT.BeamStiffness[5,5] = 9.346e+06
        XBINPUT.BeamMass[0,0] = 100
        XBINPUT.BeamMass[1,1] = 100
        XBINPUT.BeamMass[2,2] = 100
        XBINPUT.BeamMass[3,3] = 10
        XBINPUT.BeamMass[4,4] = 0.0001
        XBINPUT.BeamMass[5,5] = 0.0001
          
        "Dynamic parameters"
        XBINPUT.t0 = 0.0
        XBINPUT.tfin = 2.0
        XBINPUT.dt = 0.001
        XBINPUT.Omega = 20.0
        XBINPUT.ForceDyn[-1,2] = 60e+03
        XBINPUT.ForcingType = 'RampSin'
        XBINPUT.RampTime = 1.0
      
        XBOPTS.PrintInfo.value = False
        NonlinearDynamic.Solve_F90_steps(XBINPUT,XBOPTS)
  
        "Read from file"
        Dyn = np.loadtxt(Settings.OutputDir + Settings.OutputFileRoot + '_SOL'+\
                        str(XBOPTS.Solution.value) + '_dyn.dat')
      
          
        "Assert z-displacement at 1, 1.5 & 2 secs"
        self.assertAlmostEqual(Dyn[1000,3],DynF90[1000,3], 3,\
                                'z-displacement does not match')
        self.assertAlmostEqual(Dyn[1500,3],DynF90[1500,3], 3,\
                                'z-displacement does not match')
        self.assertAlmostEqual(Dyn[2000,3],DynF90[2000,3], 3,\
                                'z-displacement does not match')
         
     
    def test_DynTipDispRot_Py(self):
        """@brief Compare to values of Py to those of F90 steps."""
         
        "Read file from previous test (F90steps)"
        DynSteps = np.loadtxt(Settings.OutputDir+Settings.OutputFileRoot+'_F90steps'
                              +'_SOL'+ str(312) + '_dyn.dat')
         
        Settings.OutputFileRoot = 'PyBeam2'
        import NonlinearDynamic
         
        """Set up Xbopts for nonlinear static analysis defined in 
        NonlinearStatic/testcases.pdf case 1.2."""
        XBOPTS = DerivedTypes.Xbopts()
        XBOPTS.Solution.value = 312 
        XBOPTS.NumLoadSteps.value = 10
        XBOPTS.MinDelta.value = 1e-05
        XBOPTS.FollowerForce.value = False
        XBOPTS.NumGauss.value = 1
                
        """Set up Xbinput for nonlinear static analysis defined in 
        NonlinearStatic/testcases.pdf case 1.2."""
        XBINPUT = DerivedTypes.Xbinput(2,20)
        XBINPUT.BeamLength = 5.0
        XBINPUT.BeamStiffness[0,0] = 4.8e+08
        XBINPUT.BeamStiffness[1,1] = 3.231e+08
        XBINPUT.BeamStiffness[2,2] = 3.231e+08
        XBINPUT.BeamStiffness[3,3] = 1.0e+06
        XBINPUT.BeamStiffness[4,4] = 9.346e+06
        XBINPUT.BeamStiffness[5,5] = 9.346e+06
        XBINPUT.BeamMass[0,0] = 100
        XBINPUT.BeamMass[1,1] = 100
        XBINPUT.BeamMass[2,2] = 100
        XBINPUT.BeamMass[3,3] = 10
        XBINPUT.BeamMass[4,4] = 0.0001
        XBINPUT.BeamMass[5,5] = 0.0001
         
        "Dynamic parameters"
        XBINPUT.t0 = 0.0
        XBINPUT.tfin = 2.0
        XBINPUT.dt = 0.001
        XBINPUT.Omega = 20.0
        XBINPUT.ForceDyn[-1,2] = 60e+03
        XBINPUT.ForcingType = 'RampSin'
        XBINPUT.RampTime = 1.0
     
        XBOPTS.PrintInfo.value = False
        NonlinearDynamic.Solve_Py(XBINPUT,XBOPTS)
 
        "Read from file"
        Dyn = np.loadtxt(Settings.OutputDir + Settings.OutputFileRoot + '_SOL'+\
                        str(XBOPTS.Solution.value) + '_dyn.dat')
     
         
        "Assert z-displacement at 1, 1.5 & 2 secs"
        self.assertAlmostEqual(Dyn[1000,3],DynSteps[1000,3], 3,\
                                'z-displacement does not match')
        self.assertAlmostEqual(Dyn[1500,3],DynSteps[1500,3], 3,\
                                'z-displacement does not match')
        self.assertAlmostEqual(Dyn[2000,3],DynSteps[2000,3], 3,\
                                'z-displacement does not match')
        

class TestNonlinearDynamic_GolandFree(unittest.TestCase):
    """@brief Free Vibration of Goland wing."""
     
    def setUp(self):
        "Set SharPy output directory and file root"
        Settings.OutputDir = TestDir
        
    
    def tearDown(self):
        pass
    
    
    def test_DynTipDispRot_2v3noded(self):
        """@brief Check free vibration frequencies of goland wing."""

        import NonlinearDynamic # imported after clean/make process
        
        Settings.OutputFileRoot = 'PyBeamGolandFree2noded'
        
        "beam options"
        XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),\
                                     MaxIterations = ct.c_int(99),\
                                     PrintInfo = ct.c_bool(True),\
                                     NumLoadSteps = ct.c_int(25),\
                                     Solution = ct.c_int(312),\
                                     MinDelta = ct.c_double(1e-4),\
                                     NewmarkDamp = ct.c_double(1e-2))
        
        "beam inputs"
        XBINPUT = DerivedTypes.Xbinput(2,24)
        XBINPUT.BeamLength = 6.096
        XBINPUT.BeamStiffness[0,0] = 1.0e+09
        XBINPUT.BeamStiffness[1,1] = 1.0e+09
        XBINPUT.BeamStiffness[2,2] = 1.0e+09
        XBINPUT.BeamStiffness[3,3] = 0.9875e+06
        XBINPUT.BeamStiffness[4,4] = 9.77e+06
        XBINPUT.BeamStiffness[5,5] = 9.77e+08
        
        XBINPUT.BeamStiffness[:,:] = 1.0*XBINPUT.BeamStiffness[:,:]
        
        XBINPUT.BeamMass[0,0] = 35.709121
        XBINPUT.BeamMass[1,1] = 35.709121
        XBINPUT.BeamMass[2,2] = 35.709121
        XBINPUT.BeamMass[3,3] = 8.6405832
        XBINPUT.BeamMass[4,4] = 0.001
        XBINPUT.BeamMass[5,5] = 0.001
        
        "pitch-plunge coupling term"
        "b-frame coordinates"
        c = 1.8288
        
        "using skew-symmetric operator"
        cg = np.array([0.0, -0.1, 0.0])*c
        cgSkew = np.array([[   0.0, -cg[2], cg[1] ],\
                           [ cg[2],    0.0, -cg[0]],\
                           [-cg[1],  cg[0],   0.0] ])
        
        XBINPUT.BeamMass[:3,3:] = -XBINPUT.BeamMass[0,0]*cgSkew
        XBINPUT.BeamMass[3:,:3] = XBINPUT.BeamMass[:3,3:].T
        
        "set (0,5) and (5,0) to zero"
        XBINPUT.BeamMass[0,5] = 0.0
        XBINPUT.BeamMass[5,0] = 0.0
        
        
        "Dynamic parameters"
        XBINPUT.t0 = 0.0
        XBINPUT.tfin = 1.0
        XBINPUT.dt = 0.01
        XBINPUT.Omega = 0.0
        XBINPUT.ForceDyn[-1,2] = 6e03
        XBINPUT.ForcingType = '1-cos'
    
        NonlinearDynamic.Solve_Py(XBINPUT,XBOPTS)
        
        "Read from file"
        Dyn2noded = np.loadtxt(Settings.OutputDir + Settings.OutputFileRoot+'_SOL'+str(XBOPTS.Solution.value) + '_dyn.dat', \
                               skiprows=3)
        
        
        "run with 3-noded"
        Settings.OutputFileRoot = 'PyBeamGolandFree2noded'
        "beam inputs"
        XBINPUT = DerivedTypes.Xbinput(3,12)
        XBINPUT.BeamLength = 6.096
        XBINPUT.BeamStiffness[0,0] = 1.0e+09
        XBINPUT.BeamStiffness[1,1] = 1.0e+09
        XBINPUT.BeamStiffness[2,2] = 1.0e+09
        XBINPUT.BeamStiffness[3,3] = 0.9875e+06
        XBINPUT.BeamStiffness[4,4] = 9.77e+06
        XBINPUT.BeamStiffness[5,5] = 9.77e+08
        
        XBINPUT.BeamStiffness[:,:] = 1.0*XBINPUT.BeamStiffness[:,:]
        
        XBINPUT.BeamMass[0,0] = 35.709121
        XBINPUT.BeamMass[1,1] = 35.709121
        XBINPUT.BeamMass[2,2] = 35.709121
        XBINPUT.BeamMass[3,3] = 8.6405832
        XBINPUT.BeamMass[4,4] = 0.001
        XBINPUT.BeamMass[5,5] = 0.001
        
        "pitch-plunge coupling term"
        "b-frame coordinates"
        c = 1.8288
        
        "using skew-symmetric operator"
        cg = np.array([0.0, -0.1, 0.0])*c
        cgSkew = np.array([[   0.0, -cg[2], cg[1] ],\
                           [ cg[2],    0.0, -cg[0]],\
                           [-cg[1],  cg[0],   0.0] ])
        
        XBINPUT.BeamMass[:3,3:] = -XBINPUT.BeamMass[0,0]*cgSkew
        XBINPUT.BeamMass[3:,:3] = XBINPUT.BeamMass[:3,3:].T
        
        "set (0,5) and (5,0) to zero"
        XBINPUT.BeamMass[0,5] = 0.0
        XBINPUT.BeamMass[5,0] = 0.0
        
        "Dynamic parameters"
        XBINPUT.t0 = 0.0
        XBINPUT.tfin = 1.0
        XBINPUT.dt = 0.001
        XBINPUT.Omega = 0.0
        XBINPUT.ForceDyn[-1,2] = 6e03
        XBINPUT.ForcingType = '1-cos'
        
        NonlinearDynamic.Solve_Py(XBINPUT,XBOPTS)
        
        "Read from file"
        Dyn3noded = np.loadtxt(Settings.OutputDir + Settings.OutputFileRoot + '_SOL' + str(XBOPTS.Solution.value) + '_dyn.dat', \
                               skiprows=3)
        
        PlotThings = False
        if PlotThings == True:
            "plot displacements"
            "plot displacements"
            plt.figure(1)
            plt.subplot(311)
            plt.plot(Dyn3noded[:,0],Dyn3noded[:,1])
            plt.plot(Dyn2noded[:,0],Dyn2noded[:,1])
            plt.grid()
            plt.ylabel('R_1')
            #plt.ylim(4.95,5.05)
            plt.subplot(312)
            plt.plot(Dyn3noded[:,0],Dyn3noded[:,2])
            plt.plot(Dyn2noded[:,0],Dyn2noded[:,2])
            plt.grid()
            plt.ylabel('R_2')
            plt.subplot(313)
            plt.plot(Dyn3noded[:,0],Dyn3noded[:,3])
            plt.plot(Dyn2noded[:,0],Dyn2noded[:,3])
            plt.grid()
            plt.ylabel('R_3')
            plt.xlabel('time')
            plt.legend(("3-noded","2-noded")\
                       ,'lower left',\
                       shadow=True,\
                       prop={'size':12})
            plt.figure(2)
            plt.subplot(311)
            plt.plot(Dyn3noded[:,0],Dyn3noded[:,4])
            plt.plot(Dyn2noded[:,0],Dyn2noded[:,4])
            plt.grid()
            plt.ylabel('Psi_1')
            #plt.ylim(4.95,5.05)
            plt.subplot(312)
            plt.plot(Dyn3noded[:,0],Dyn3noded[:,5])
            plt.plot(Dyn2noded[:,0],Dyn2noded[:,5])
            plt.grid()
            plt.ylabel('Psi_2')
            plt.subplot(313)
            plt.plot(Dyn3noded[:,0],Dyn3noded[:,6])
            plt.plot(Dyn2noded[:,0],Dyn2noded[:,6])
            plt.grid()
            plt.ylabel('Psi_3')
            plt.xlabel('time')
            plt.legend(("3-noded","2-noded")\
                       ,'lower left',\
                       shadow=True,\
                       prop={'size':12})
            plt.show()
        
        for iTime in range(Dyn3noded.shape[0]):
            "Assert z-displacement"
            Delz = 0.001*np.max(Dyn3noded[:,3])
            self.assertAlmostEqual(Dyn2noded[iTime,3],Dyn3noded[iTime,3],\
                                   places = None,\
                                   msg = 'z-displacement does not match',\
                                   delta=Delz)
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    #unittest.main()
    suite1 = unittest.TestLoader().loadTestsFromTestCase(\
                TestNonlinearDynamic_v_Executable)
    suite2 = unittest.TestLoader().loadTestsFromTestCase(\
                TestNonlinearDynamic_v_PalaciosCesnik)
    suite3 = unittest.TestLoader().loadTestsFromTestCase(\
                TestNonlinearDynamic_GolandFree)
    alltests = unittest.TestSuite([suite1, suite2, suite3])
    #alltests.run(unittest.defaultTestResult())
    TestRunner = unittest.TextTestRunner(verbosity=2) # creates a test runner
    TestRunner.run(suite2) #run a single suite
    #TestRunner.run(alltests) #run all tests