'''
Created on 14 Jul 2015

@author: rjs10
'''
import unittest
import SharPySettings as Settings
import numpy as np
import ctypes as ct
from PyAero.UVLM.Utils.Linear import genSSuvlm
from PyAero.UVLM.Utils.UVLMLib import Cpp_KJForces, Cpp_Solver_VLM, Cpp_KJForces_vC
from PyAero.UVLM.Utils.DerivedTypesAero import VMopts
from Misc import isodd

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
        m=2
        n=2
        mW=3
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
         
        # initialise nonlinear options
        VMOPTS = VMopts(m, n, False, # image methods
                        mW,
                        True, # steady
                        True, # KJ is true
                        True,
                        delS,
                        False, 
                        1) #numCores
         
        # solve nonlinear
        Cpp_Solver_VLM(zeta0, zetaPri0, nu, zetaW0, VMOPTS, f0, gam0, gamW0)
         
        # check lift curve slope
        self.assertAlmostEqual((sum((f0[2::3])/10000.0))/aoa,2*np.pi,1)
         
        # test independent force calculation
        f0i = np.zeros((3*len(chords)*len(spans)))
        gamPri0=np.zeros((m*n))
        gam_tm1=gam0-gamPri0
        Cpp_KJForces(zeta0, gam0, zetaW0, gamW0, zetaPri0, nu, VMOPTS, gam_tm1, f0i)
        self.assertAlmostEqual((sum((f0i[2::3])/10000.0))/aoa,2*np.pi,1)
        self.assertTrue(sum(f0i[2::3]) == sum(f0[2::3]))      
         
    def test_linVnln(self):
        """Test ss UVLM output equations against full nonlinear calcs."""
          
        # Init steady problem at 1 deg AoA
        aoa = 1*np.pi/180.0
        V = 1
        m=2
        n=2
        mW=3
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
        gamPri0=np.zeros((m*n))
        nu = np.zeros((3*len(chords)*len(spans))) # atmospheric velocity
        nu[0::3]=V
         
        # to be solved 
        gam0=np.zeros((m*n))
        gamW0=np.zeros((mW*n))
        f0 = np.zeros((3*len(chords)*len(spans))) # forces
         
        # initialise nonlinear options
        VMOPTS = VMopts(m, n, False, # image methods
                        mW,
                        True, # steady
                        True, # KJ is true
                        True,
                        delS,
                        False, 
                        1) #numCores
         
        # solve nonlinear
        Cpp_Solver_VLM(zeta0, zetaPri0, nu, zetaW0, VMOPTS, f0, gam0, gamW0)
         
         
        # generate linear output eqs at x0, u0
        foo1, foo2, foo3, C, D = genSSuvlm(gam0, gamW0, gamPri0, zeta0, zetaW0, zetaPri0, nu, m, n, mW, delS)
        del foo1, foo2, foo3
         
        # init delta vectors for testing
        dX = np.zeros((2*m*n+mW*n))
        dU = np.zeros((9*len(chords)*len(spans)))
         
        # variations in gamma--------------------------------------------------
        for i in range(m*n):
            dX[i] = 0.00002/(m*n)*(i+1)
        gamPdX = gam0+dX[0:m*n]
        f = np.zeros_like(f0)
        gam_tm1=gam0-delS*gamPri0
        Cpp_KJForces(zeta0, gamPdX, zetaW0, gamW0, zetaPri0, nu, VMOPTS, gam_tm1, f)
        f_dGamma = f - f0
        dfApprox = np.dot(C,dX)+np.dot(D,dU)
         
        # lift
        dLexact = sum(f_dGamma[2::3])
        dLapprox = sum(dfApprox[2::3])
        self.assertLess(np.abs((dLapprox-dLexact)/dLexact), 0.001)
         
        # drag
        dDexact = sum(f_dGamma[0::3])
        dDapprox = sum(dfApprox[0::3])
        self.assertLess(np.abs((dDapprox-dDexact)/dDexact), 0.001)
         
        # side force
        dSexact = sum(f_dGamma[1::3])
        dSapprox = sum(dfApprox[1::3])
        self.assertLess(dSexact, 1e-7)
        self.assertLess(dSapprox, 1e-7)
        
        # variations in gamma w ------------------------------------------------
        f[:]=0.0
        dX[:]=0.0
        for i in range(m*n,mW*n):
            dX[i] = 0.00002/(mW*n)*(i+1)
        gamWpDx = gamW0 + dX[m*n:m*n+mW*n]
        Cpp_KJForces(zeta0, gam0, zetaW0, gamWpDx, zetaPri0, nu, VMOPTS, gam_tm1, f)
        f_dGammaW = f - f0
        dfApprox = np.dot(C,dX)+np.dot(D,dU)
         
        # lift
        dLexact = sum(f_dGammaW[2::3])
        dLapprox = sum(dfApprox[2::3])
        self.assertLess(np.abs((dLapprox-dLexact)/dLexact), 0.001)
         
        # drag
        dDexact = sum(f_dGammaW[0::3])
        dDapprox = sum(dfApprox[0::3])
        self.assertLess(np.abs((dDapprox-dDexact)/dDexact), 0.001)
         
        # side force
        dSexact = sum(f_dGammaW[1::3])
        dSapprox = sum(dfApprox[1::3])
        self.assertLess(dSexact, 1e-7)
        self.assertLess(dSapprox, 1e-7)
         
        # variations in gamPri ------------------------------------------------
        f[:]=0.0
        dX[:]=0.0
        for i in range(m*n+mW*n,2*m*n+mW*n):
            dX[i] = 0.002/(m*n)*(i+1)
        gam_tm1pDx = gam0-2.0*delS*dX[m*n+mW*n:]
        VMOPTS.Steady = ct.c_bool(False)
        Cpp_KJForces(zeta0, gam0, zetaW0, gamW0, zetaPri0, nu, VMOPTS, gam_tm1pDx, f)
        f_dGamPri = f - f0
        dfApprox = np.dot(C,dX)+np.dot(D,dU)
         
        # lift
        dLexact = sum(f_dGamPri[2::3])
        dLapprox = sum(dfApprox[2::3])
        self.assertLess(np.abs((dLapprox-dLexact)/dLexact), 0.001)
         
        # drag
        dDexact = sum(f_dGamPri[0::3])
        dDapprox = sum(dfApprox[0::3])
        self.assertLess(np.abs((dDapprox-dDexact)/dDexact), 0.001)
         
        # side force
        dSexact = sum(f_dGamPri[1::3])
        dSapprox = sum(dfApprox[1::3])
        self.assertLess(dSexact, 1e-7)
        self.assertLess(dSapprox, 1e-7)
         
        # variations in zetaPri ------------------------------------------------
        f[:]=0.0
        dX[:]=0.0
        for i in range(3*(m+1)*(n+1)):
            dU[i] = 0.002/(3*(m+1)*(n+1)*(i+1))
        zetaPriPdX = zetaPri0+dU[0:3*(m+1)*(n+1)]
        VMOPTS.Steady = ct.c_bool(True)
        Cpp_KJForces(zeta0, gam0, zetaW0, gamW0, zetaPriPdX, nu, VMOPTS, gam_tm1, f)
        f_dZetaPri = f - f0
        dfApprox = np.dot(C,dX)+np.dot(D,dU)
         
        # lift
        dLexact = sum(f_dZetaPri[2::3])
        dLapprox = sum(dfApprox[2::3])
        self.assertLess(np.abs((dLapprox-dLexact)/dLexact), 0.001)
         
        # drag
        dDexact = sum(f_dZetaPri[0::3])
        dDapprox = sum(dfApprox[0::3])
        self.assertLess(np.abs((dDapprox-dDexact)/dDexact), 0.001)
         
        # side force
        dSexact = sum(f_dZetaPri[1::3])
        dSapprox = sum(dfApprox[1::3])
        self.assertLess(np.abs((dSapprox-dSexact)/dSexact), 0.001)
         
        # variations in nu -----------------------------------------------------
        f[:]=0.0
        dU[:]=0.0
        for i in range(6*(m+1)*(n+1),9*(m+1)*(n+1)):
            dU[i] = 0.002/(3*(m+1)*(n+1)*(i+1))
        nuPdX = nu + dU[6*(m+1)*(n+1):]
        Cpp_KJForces(zeta0, gam0, zetaW0, gamW0, zetaPri0, nuPdX, VMOPTS, gam_tm1, f)
        f_dNu = f - f0
        dfApprox = np.dot(C,dX)+np.dot(D,dU)
         
        # lift
        dLexact = sum(f_dNu[2::3])
        dLapprox = sum(dfApprox[2::3])
        self.assertLess(np.abs((dLapprox-dLexact)/dLexact), 0.001)
         
        # drag
        dDexact = sum(f_dNu[0::3])
        dDapprox = sum(dfApprox[0::3])
        self.assertLess(np.abs((dDapprox-dDexact)/dDexact), 0.001)
         
        # side force
        dSexact = sum(f_dNu[1::3])
        dSapprox = sum(dfApprox[1::3])
        self.assertLess(np.abs((dSapprox-dSexact)/dSexact), 0.001)
        

    def test_linVnln_dZeta(self):
        """Test ss UVLM output equations against full nonlinear calcs."""
         
        # Init steady problem at 1 deg AoA
        aoa = 1.0*np.pi/180.0
        V = 1
        m=4
        n=3
        mW=11
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
        gamPri0=np.zeros((m*n))
        nu = np.zeros((3*len(chords)*len(spans))) # atmospheric velocity
        nu[0::3]=V
        
        # to be solved 
        gam0=np.zeros((m*n))
        gamW0=np.zeros((mW*n))
        f0 = np.zeros((3*len(chords)*len(spans))) # forces
        
        # initialise nonlinear options
        VMOPTS = VMopts(m, n, False, # image methods
                        mW,
                        True, # steady
                        True, # KJ is true
                        True,
                        delS,
                        False, 
                        1) #numCores
        
        # solve nonlinear
        Cpp_Solver_VLM(zeta0, zetaPri0, nu, zetaW0, VMOPTS, f0, gam0, gamW0)
        
        # set to unsteady mode and calculate with zetaPri0
        f0[:]=0.0
        VMOPTS.Steady = ct.c_bool(False)
        gamPri0=0.01*np.ones_like(gamPri0)
        gam_tm1=gam0-2.0*delS*gamPri0
        Cpp_KJForces(zeta0, gam0, zetaW0, gamW0, zetaPri0, nu, VMOPTS, gam_tm1, f0)
        
        # generate linear output eqs at x0, u0
        foo1, foo2, foo3, C, D = genSSuvlm(gam0, gamW0, gamPri0, zeta0, zetaW0, zetaPri0, nu, m, n, mW, delS)
        del foo1, foo2, foo3
        
        # init delta vectors for testing
        dX = np.zeros((2*m*n+mW*n))
        dU = np.zeros((9*len(chords)*len(spans)))
        
        # variations in zeta -----------------------------------------------------
        f = np.zeros_like(f0)
        for i in range(3*(m+1)*(n+1),6*(m+1)*(n+1)):
            if isodd(i):
                dU[i] = 0.1/(m*3*(m+1)*(n+1)*(i+1))
            else:
                dU[i] = -0.1/(m*3*(m+1)*(n+1)*(i+1))
        zetaPdX = zeta0 + dU[3*(m+1)*(n+1):6*(m+1)*(n+1)]
        gam_tm1=gam0-2.0*delS*gamPri0
        Cpp_KJForces(zetaPdX, gam0, zetaW0, gamW0, zetaPri0, nu, VMOPTS, gam_tm1, f)
        
        # calculate diffs
        f_dZeta = f - f0
        dfApprox = np.dot(C,dX)+np.dot(D,dU)
        
        # lift
        dLexact = sum(f_dZeta[2::3])
        dLapprox = sum(dfApprox[2::3])
        self.assertLess(np.abs((dLapprox-dLexact)/dLexact), 0.01)
        
        # drag
        dDexact = sum(f_dZeta[0::3])
        dDapprox = sum(dfApprox[0::3])
        self.assertLess(np.abs((dDapprox-dDexact)/dDexact), 0.01)
        
        # side force
        dSexact = sum(f_dZeta[1::3])
        dSapprox = sum(dfApprox[1::3])
        self.assertLess(np.abs((dSapprox-dSexact)/dSexact), 0.01)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()