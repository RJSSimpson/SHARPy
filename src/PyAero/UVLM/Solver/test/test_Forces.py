'''
Created on 14 Jul 2015

@author: rjs10
'''
import unittest
import SharPySettings as Settings
import numpy as np
from PyAero.UVLM.Utils.Linear import genSSuvlm
from PyAero.UVLM.Utils.UVLMLib import Cpp_KJForces, Cpp_Solver_VLM, Cpp_KJForces_vC
from PyAero.UVLM.Utils.DerivedTypesAero import VMopts

TestDir = (Settings.SharPyProjectDir + 'output/tests/PyAero/UVLM/' 
           + 'Forces/')

class Test_linearForces(unittest.TestCase):


    def setUp(self):
        # Set SharPy output directory and file root.
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'linVnln'


    def tearDown(self):
        pass
    
    def test_nlnVsolver(self):
        """Test nln force calcs against full solver solution on aerofoil
        problem."""
        # Init steady problem at 1 deg AoA
        aoa = 1*np.pi/180.0
        V = 1
        m=1
        n=1
        mW=1
        delS=1
        chords = np.linspace(-1.0, 0.0, m+1, True)
        chordsW = np.linspace(0.0, 10000.0, mW+1, True)
        spans = np.linspace(-10000, 10000, n+1, True)
        zeta0=np.zeros(3*len(chords)*len(spans))
        zetaW0=np.zeros(3*len(chordsW)*len(spans))
        kk=0
        for c in chords:
            for s in spans:
                zeta0[3*kk]=np.cos(aoa)*c
                zeta0[3*kk+1]=s
                zeta0[3*kk+2]=-np.sin(aoa)*c
                kk=kk+1
            # end for s
        # end for c
        kk=0
        for c in chordsW:
            for s in spans:
                zetaW0[3*kk]=c
                zetaW0[3*kk+1]=s
                kk=kk+1
            # end for s
        # end for c
        zetaPri0 = np.zeros((3*len(chords)*len(spans)))
        nu = np.zeros((3*len(chords)*len(spans))) # atmospheric velocity
        nu[0::3]=V
        
        # to be solved 
        gam0=np.zeros((m*n))
        gamW0=np.zeros((mW*n))
        f0 = np.zeros((3*len(chords)*len(spans))) # forces
        
        VMOPTS = VMopts(m, n, False, # image methods
                        mW,
                        True, # steady
                        True, # KJ is true
                        True,
                        delS,
                        False, 
                        1) #numCores
        
        Cpp_Solver_VLM(zeta0, zetaPri0, nu, zetaW0, VMOPTS, f0, gam0, gamW0)
        self.assertAlmostEqual((sum((f0[2::3])/10000.0))/aoa,2*np.pi,1)
        print("\n steady D = ", sum(f0[0:3])/10000.0)
        
        # test independent force calculation
        f0i = np.zeros((3*len(chords)*len(spans)))
        gamPri0=np.zeros((m*n))
        gam_tm1=gam0-gamPri0
        Cpp_KJForces_vC(zeta0, gam0, zetaW0, gamW0, zetaPri0, nu, VMOPTS, gam_tm1, f0i)
        self.assertAlmostEqual((sum((f0[2::3])/10000.0))/aoa,2*np.pi,1)
        print("\n steady D = ", sum(f0i[0:3])/10000.0)
        
        # generate linear output eqs at x0, u0
        foo1, foo2, foo3, C, D = genSSuvlm(gam0, gamW0, gamPri0, zeta0, zetaW0, zetaPri0, nu, m, n, mW, delS)
        del foo1, foo2, foo3
        
        # innit delta vectors for testing
        dX = np.zeros((2*m*n+mW*n))
        dU = np.zeros((9*len(chords)*len(spans)))
        
        # variations in gamma
        for i in range(m*n):
            dX[i] = 0.00002/(m*n)*(i+1)
        gamPdX = gam0+dX[0:m*n]
        f0idGamma = np.zeros_like(f0i)
        Cpp_KJForces_vC(zeta0, gamPdX, zetaW0, gamW0, zetaPri0, nu, VMOPTS, gam_tm1, f0idGamma)
        f_dGamma = f0idGamma - f0i
        dfApprox = np.dot(C,dX)+np.dot(D,dU)
        print("\n df_dGamma:\n",f_dGamma)
        print("\n dfApprox:\n",dfApprox)
        print("\n diff:\n",f_dGamma-dfApprox)
        
        dLexact = sum(f_dGamma[2::3])
        dLapprox = sum(dfApprox[2::3])
        print("\n dL (exact) = ",dLexact,"\t dL (approx) = ", dLapprox)
        
        dDexact = sum(f_dGamma[0::3])
        dDapprox = sum(dfApprox[0::3])
        print("\n dD (exact) = ",dDexact,"\t dD (approx) = ", dDapprox)
        
        dSexact = sum(f_dGamma[1::3])
        dSapprox = sum(dfApprox[1::3])
        print("\n dS (exact) = ",dSexact,"\t dS (approx) = ", dSapprox)
        
        print("\n\n---------------------------------------------\n C (gamma):\n",
              C[:,m*n])
       
        # variations in gamma w
#         dX[:]=0.0
#         dX[m*n:m*n+mW*n] = np.random.random((mW*n))/10000.0
#         gamWpDx = gamW0 + dX[m*n:m*n+mW*n]
#         f0idGammaW = np.zeros_like(f0i)
#         Cpp_KJForces(zeta0, gam0, zetaW0, gamWpDx, zetaPri0, nu, VMOPTS, gam_tm1, f0idGammaW)
#         f_dGammaW = f0idGammaW - f0i
#         print(f_dGammaW)
#         print(np.dot(C,dX)+np.dot(D,dU))
#         print("\n diff:\n",f_dGammaW-np.dot(C,dX)-np.dot(D,dU))
        
        
#     def test_linVnln(self):
#         """Test ss UVLM output equations against full nonlinear calcs."""
#         
#         # define test case
#         m=1
#         n=1
#         mW=1
#         delS=1
#         gam0=np.ones((m*n))
#         gamW0=np.ones((mW*n))
#         gamPri0=0.5*np.ones((m*n))
#         chords = np.linspace(0.0, 1.0, m+1, True)
#         chordsW = np.linspace(1.0, 3.0, mW+1, True)
#         spans = np.linspace(0.0, 1.0, n+1, True)
#         zeta0=np.zeros(3*len(chords)*len(spans))
#         zetaW0=np.zeros(3*len(chordsW)*len(spans))
#         zetaPri0 = np.ones((3*len(chords)*len(spans)))
#         nu0 = np.ones((3*len(chords)*len(spans)))
#         kk=0
#         for c in chords:
#             for s in spans:
#                 zeta0[3*kk]=c
#                 zeta0[3*kk+1]=s
#                 kk=kk+1
#             # end for s
#         # end for c
#         kk=0
#         for c in chordsW:
#             for s in spans:
#                 zetaW0[3*kk]=c
#                 zetaW0[3*kk+1]=s
#                 kk=kk+1
#             # end for s
#         # end for c
#         
#         # generate linear model
#         C = genSSuvlm(gam0, gamW0, gamPri0, zeta0, zetaW0, zetaPri0, nu0, m, n, mW, delS)[3]
#         D = genSSuvlm(gam0, gamW0, gamPri0, zeta0, zetaW0, zetaPri0, nu0, m, n, mW, delS)[4]
#         
#         # get reference forces
#         f0 = np.zeros((3*len(chords)*len(spans))) # forces
#         VMOPTS = VMopts(m, n, False, # image methods
#                         mW,
#                         False,
#                         True,
#                         False,
#                         delS,
#                         False, 
#                         1) #numCores
#         
#         # get gamma at previous time step
#         gam_tm1=gam0-delS*gamPri0
#         
#         # claculate forces
#         Cpp_KJForces(zeta0, gam0, zetaW0, gamW0, zetaPri0, nu0, VMOPTS, gam_tm1, f0)
#         #print(f0)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()