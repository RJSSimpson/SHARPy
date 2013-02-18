'''
Created on 15 Feb 2013

@author: rjs10
'''
import unittest
import SharPySettings as Settings
import numpy as np
import ctypes as ct
import DerivedTypes
import DerivedTypesAero
from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps
import PyCoupled.Coupled_NlnStatic as Solver
from _pyio import open


TestDir = Settings.SharPyProjectDir + 'SharPy/src/PyCoupled/' \
           + 'test/NlnStatic/'

class Test_CantHALE_v_Murua2012_Smith2001(unittest.TestCase):


    def setUp(self):
        "Set SharPy output directory and file root"
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'CantHALE4'


    def tearDown(self):
        pass


    def test_BendingDisp(self):
        """@brief Test tip disp and plot bening against Murua and Smith."""
        
        "beam options"
        XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),\
                                     MaxIterations = ct.c_int(50),\
                                     PrintInfo = ct.c_bool(True),\
                                     NumLoadSteps = ct.c_int(25),\
                                     Solution = ct.c_int(112),\
                                     MinDelta = ct.c_double(1e-4))
        
        "beam inputs"
        XBINPUT = DerivedTypes.Xbinput(3,20)
        XBINPUT.BeamLength = 16.0
        XBINPUT.BeamStiffness[0,0] = 1.0e+09
        XBINPUT.BeamStiffness[1,1] = 1.0e+09
        XBINPUT.BeamStiffness[2,2] = 1.0e+09
        XBINPUT.BeamStiffness[3,3] = 1.0e+04
        XBINPUT.BeamStiffness[4,4] = 2.0e+04
        XBINPUT.BeamStiffness[5,5] = 5.0e+06
        
        XBINPUT.BeamStiffness[:,:] = 1.0*XBINPUT.BeamStiffness[:,:]
        
        XBINPUT.BeamMass[0,0] = 0.75
        XBINPUT.BeamMass[1,1] = 0.75
        XBINPUT.BeamMass[2,2] = 0.75
        XBINPUT.BeamMass[3,3] = 0.1
        XBINPUT.BeamMass[4,4] = 0.001
        XBINPUT.BeamMass[5,5] = 0.001
        
        
        "aero options"
        VMOPTS = DerivedTypesAero.VMopts(M = 1,\
                                         N = XBINPUT.NumNodesTot - 1 ,\
                                         ImageMethod = True,\
                                         Steady = True,\
                                         KJMeth = False)
        
        "aero inputs"
        VMINPUT = DerivedTypesAero.VMinput(c = 1.0, b = XBINPUT.BeamLength,\
                                           U_mag = 25.0,\
                                           alpha = 4.0*np.pi/180.0,\
                                           theta = 0.0)
        
        "aeroelastic opts"
        "Density and gravity due to US standard atmosphere at 20km"
        AELAOPTS = AeroelasticOps(0.0,0.0,0.08891)
        
               
        PosDefor, PsiDefor = Solver.Solve_Py(XBINPUT,XBOPTS,\
                                             VMOPTS,VMINPUT,\
                                             AELAOPTS)
        
        self.assertAlmostEqual(PosDefor[-1,2], 5.336874, None, \
                               'Tip deflection wrong.', 0.02)
        
        convergence = False
        if convergence == True:
            
            VMOPTS.KJMeth = ct.c_bool(True)
            
            FileMesh = open(TestDir + 'ConvergenceMesh.dat', "w")
            FileKJ = open(TestDir + 'ConvergenceKJ.dat', "w")
            FileKP = open(TestDir + 'ConvergenceKP.dat', "w")
            
            NumElemsArr = [1,5,10,15,20,30]
            MSectArr = [1,2,5,10,15]
            
            "write mesh data corresponding to test"
            FileMesh.write("Num Nodes per Elem = %d\n" % (XBINPUT.NumNodesElem))
            FileMesh.write("Num Nodes Array = [")
            for Bla in NumElemsArr:
                FileMesh.write("%d  " %(Bla))
            FileMesh.write("]\n")
            FileMesh.write("Chordwise Panels Array = [")
            for Bla in MSectArr:
                FileMesh.write("%d  " %(Bla))
            FileMesh.write("]\n")
            FileMesh.close()
            
            for NumElems in NumElemsArr:
                for MSect in MSectArr:
                    
                    "beam options"
                    XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),\
                                                 MaxIterations = ct.c_int(50),\
                                                 PrintInfo = ct.c_bool(True),\
                                                 NumLoadSteps = ct.c_int(25),\
                                                 Solution = ct.c_int(112),\
                                                 MinDelta = ct.c_double(1e-4))
                    
                    "beam inputs"
                    XBINPUT = DerivedTypes.Xbinput(3,NumElems)
                    XBINPUT.BeamLength = 16.0
                    XBINPUT.BeamStiffness[0,0] = 1.0e+09
                    XBINPUT.BeamStiffness[1,1] = 1.0e+09
                    XBINPUT.BeamStiffness[2,2] = 1.0e+09
                    XBINPUT.BeamStiffness[3,3] = 1.0e+04
                    XBINPUT.BeamStiffness[4,4] = 2.0e+04
                    XBINPUT.BeamStiffness[5,5] = 5.0e+06
                    
                    XBINPUT.BeamStiffness[:,:] = 1.0*XBINPUT.BeamStiffness[:,:]
                    
                    XBINPUT.BeamMass[0,0] = 0.75
                    XBINPUT.BeamMass[1,1] = 0.75
                    XBINPUT.BeamMass[2,2] = 0.75
                    XBINPUT.BeamMass[3,3] = 0.1
                    XBINPUT.BeamMass[4,4] = 0.001
                    XBINPUT.BeamMass[5,5] = 0.001
                    
                    
                    "aero options"
                    VMOPTS = DerivedTypesAero.VMopts(M = MSect,\
                                                     N = XBINPUT.NumNodesTot - 1 ,\
                                                     ImageMethod = True,\
                                                     Steady = True,\
                                                     KJMeth = True)
                    
                    "aero inputs"
                    VMINPUT = DerivedTypesAero.VMinput(c = 1.0, b = XBINPUT.BeamLength,\
                                                       U_mag = 25.0,\
                                                       alpha = 4.0*np.pi/180.0,\
                                                       theta = 0.0)
                    
                    "aeroelastic opts"
                    "Density and gravity due to US standard atmosphere at 20km"
                    AELAOPTS = AeroelasticOps(0.0,0.0,0.08891)
                    
                           
                    PosDefor, PsiDefor = Solver.Solve_Py(XBINPUT,XBOPTS,\
                                                         VMOPTS,VMINPUT,\
                                                         AELAOPTS)
                    
                    "write tip deflection to file"
                    FileKJ.write("%12.5e  " %(PosDefor[-1,2]))
                    if MSect == MSectArr[-1]:
                        FileKJ.write("\n")
                    
                    
                    "do K&P"
                    VMOPTS.KJMeth = ct.c_bool(False)
                    
                    PosDefor, PsiDefor = Solver.Solve_Py(XBINPUT,XBOPTS,\
                                                         VMOPTS,VMINPUT,\
                                                         AELAOPTS)
                    
                    "write tip deflection to file"
                    FileKP.write("%12.5e  "%(PosDefor[-1,2]))
                    if MSect == MSectArr[-1]:
                        FileKP.write("\n")
                    
                # END for MSect
            # END for NumNodes
            FileKJ.close()
            FileKP.close()
        # END if convergence
                        
                    
                    
                    


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()