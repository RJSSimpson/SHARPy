'''
Created on 25 Feb 2013

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
from PyCoupled.Coupled_NlnDynamic import VMCoupledUnstInput, Solve_Py

TestDir = Settings.SharPyProjectDir + 'SharPy/src/PyCoupled/' \
           + 'test/NlnDynamic/'

Settings.PlotTec = True

class Test_GolandFlutter(unittest.TestCase):


    def setUp(self):
        "Set SharPy output directory"
        Settings.OutputDir = TestDir + 'GolandFlutter/Convergence/'


    def tearDown(self):
        pass


    def test_U_160_165_170(self):
        
        
#        NumNodesElemArr = [3, 2]
        NumNodesElemArr = [2]
        
#        NumElemsArr = [12,24,36]
        NumElemsArr = [12,24,36]

        #iMArr = [6, 12, 18]
        iMArr = [6]
        
#       U = [160.0, 161.0, 162.0, 163.0, 164.0, 165.0, 166.0, 167.0, 168.0, 169.0, \
#             170.0]
        U = [164.0]
        
        Tights = [False]
        Rollup = False
        KJMeth = True
        ImpStart = False
        
        
        for NumNodesElem in NumNodesElemArr:
            for NumElems in NumElemsArr:
                for iM in iMArr:
                    for Tight in Tights:
                        
                        "beam options"
                        XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),\
                                                     MaxIterations = ct.c_int(50),\
                                                     PrintInfo = ct.c_bool(True),\
                                                     NumLoadSteps = ct.c_int(25),\
                                                     Solution = ct.c_int(312),\
                                                     MinDelta = ct.c_double(1e-5),\
                                                     NewmarkDamp = ct.c_double(1.0e-4))
                        
                        if NumNodesElem == 2:
                            NumElemsHere = NumElems
                        elif NumNodesElem == 3:
                            NumElemsHere = int(0.5*NumElems)
                        else:
                            print("How many nodes per elem!!?")
                            return
                        
                        "beam inputs"
                        XBINPUT = DerivedTypes.Xbinput(NumNodesElem,NumElemsHere)
                        XBINPUT.BeamLength = 6.096
                        XBINPUT.BeamStiffness[0,0] = 1.0e+09
                        XBINPUT.BeamStiffness[1,1] = 1.0e+09
                        XBINPUT.BeamStiffness[2,2] = 1.0e+09
                        XBINPUT.BeamStiffness[3,3] = 0.9875e+06
                        XBINPUT.BeamStiffness[4,4] = 9.77e+06
                        XBINPUT.BeamStiffness[5,5] = 9.77e+08
                        
                        XBINPUT.BeamStiffness[:,:] = 1.0*XBINPUT.BeamStiffness[:,:]
                        
                        XBINPUT.BeamMass[0,0] = 35.71
                        XBINPUT.BeamMass[1,1] = 35.71
                        XBINPUT.BeamMass[2,2] = 35.71
                        XBINPUT.BeamMass[3,3] = 8.64
                        XBINPUT.BeamMass[4,4] = 0.0001
                        XBINPUT.BeamMass[5,5] = 0.0001
                        
                        "pitch-plunge coupling term"
                        "b-frame coordinates"
                        c = 1.8288
                        
                        "EA position in Theo sectional coords"
                        ElasticAxis = -0.34
                        
                        "or using skew-symmetric operator"
                        cg = np.array([0.0, -0.1, 0.0])*c
                        cgSkew = np.array([[   0.0, -cg[2], cg[1] ],\
                                           [ cg[2],    0.0, -cg[0]],\
                                           [-cg[1],  cg[0],   0.0] ])
                        
                        XBINPUT.BeamMass[:3,3:] = -XBINPUT.BeamMass[0,0]*cgSkew
                        XBINPUT.BeamMass[3:,:3] = XBINPUT.BeamMass[:3,3:].T
                        
                        "set (0,5) and (5,0) to zero"
                        XBINPUT.BeamMass[0,5] = 0.0
                        XBINPUT.BeamMass[5,0] = 0.0
                        
                        "unsteady parameters"
                        OmegaF = 70.0 #approx flutter freq
                        
                        
                        "loop through U"
                        for iU in range(len(U)):
                            U_mag = U[iU]
                            k = OmegaF*c/(2*U_mag) #get expected reduced freq
                            M = iM#int(50.0*k/np.pi) + 1#iM #get adequate panelling
                            DelTime = c/(U_mag*M) #get resulting DelTime
                            
                            "beam"
                            XBINPUT.dt = DelTime
                            XBINPUT.t0 = 0.0
                            XBINPUT.tfin = 1.0
                            
                            "aero"
                            WakeLength = 15.0*c
                            Mstar = int(WakeLength/(DelTime*U_mag))
                            
                            
                            "aero options"
                            VMOPTS = DerivedTypesAero.VMopts(M = M,\
                                                             N = XBINPUT.NumNodesTot - 1 ,\
                                                             ImageMethod = True,\
                                                             Mstar = Mstar,\
                                                             Steady = False,\
                                                             KJMeth = KJMeth,\
                                                             NewAIC = True,\
                                                             DelTime = DelTime,\
                                                             Rollup = Rollup,\
                                                             NumCores = 4)
                            
                            "aero inputs"
                            VMINPUT = DerivedTypesAero.VMinput(c = c, b = XBINPUT.BeamLength,\
                                                               U_mag = U_mag,\
                                                               alpha = 0.05*np.pi/180.0,\
                                                               theta = 0.0,\
                                                               ZetaDotTest = 0.0,\
                                                               WakeLength = WakeLength)
                            
                            "unsteady aero inputs"
                            VelA_G = np.array(([0.0,0.0,0.0]))
                            OmegaA_G = np.array(([0.0,0.0,0.0]))
                            VMUNST = VMCoupledUnstInput(VMOPTS, VMINPUT, WakeLength, 123.4, WakeLength/c,\
                                                      VelA_G, OmegaA_G)
                            
                             
                            
                            "aeroelastic opts"
                            "Density and gravity"
                            AELAOPTS = AeroelasticOps(ElasticAxis,0.0,1.02,\
                                                      Tight = Tight,\
                                                      ImpStart = ImpStart)
                            
                            
                            "set file root"
                            Settings.OutputFileRoot = 'Goland' + \
                                                      '_NumNE' + str(NumNodesElem) + \
                                                      '_NumEl' + str(NumElemsHere) + \
                                                      '_M' + str(iM) + \
                                                      '_U' + str(U_mag)
                                                      
                            if AELAOPTS.Tight == True:
                                Settings.OutputFileRoot += '_Tight'
                                
                            if Rollup == True:
                                Settings.OutputFileRoot += '_Rollup'
                                
                            if KJMeth == True:
                                Settings.OutputFileRoot += '_KJMeth'
                                
                            if ImpStart == True:
                                Settings.OutputFileRoot += '_ImpStart'
                            
                            "solve"
                            Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,VMUNST,AELAOPTS)
                    
                    # END for iU
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()