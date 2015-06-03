'''
Created on 1 Jun 2015

@author: rjs10
'''
import unittest
import SharPySettings as Settings
import numpy as np
import ctypes as ct
from UVLMLib import Cpp_AIC, Cpp_dAgamma0_dZeta
from scipy.io import savemat

TestDir = (Settings.SharPyProjectDir + 'output/tests/PyAero/UVLM/' 
           + 'Matrices/')

class Test_AIC(unittest.TestCase):
    """@brief Test AIC matrix elements."""

    def setUp(self):
        # Set SharPy output directory and file root.
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'AIC'


    def tearDown(self):
        pass


    def test_AIC_unit(self):
        """Compare AIC matrix to trusted results."""
        #Init AIC
        m=1
        n=1
        k=m*n
        AIC = np.zeros((k,k))
        # initialize grid for 3D solver (Unit VR)
        chords = (0 , 1)
        spans = (0, 1)
        zeta=np.zeros(3*len(chords)*len(spans))
        k=0
        for c in chords:
            for s in spans:
                zeta[3*k]=c
                zeta[3*k+1]=s
                k=k+1
        # call AIC matrix function
        Cpp_AIC(zeta, 1, 1, zeta, 1, 1, AIC)
        self.assertAlmostEqual(AIC[0,0],-0.900316,6)
    
    def test_AIC_unit2DM10(self):
        """Compare AIC matrix to VLM.cpp AIC matrix"""
        #Init AIC
        m=10
        n=1
        k=m*n
        AIC = np.zeros((k,k))
        # initialize grid for 3D solver (Unit VR)
        chords = np.linspace(0.0, 1.0, m+1, True)
        spans = (-1000, 1000)
        zeta=np.zeros(3*len(chords)*len(spans))
        kk=0
        for c in chords:
            for s in spans:
                zeta[3*kk]=c
                zeta[3*kk+1]=s
                kk=kk+1
        # call AIC matrix function
        Cpp_AIC(zeta, m, n, zeta, m, n, AIC)
        # load data for comparison
        data = np.loadtxt(TestDir + "unit2DM10.dat")
        # check all entries
        for i in range(k-1):
            for j in range(k-1):
                self.assertAlmostEqual(AIC[i,j],data[i,j],4)
            # end for j
        # end for i
        
class dAgamma0_dzeta(unittest.TestCase):
    """@brief Test non-zero reference gamma matrix, dAgamma0_dzeta."""
    
    def setUp(self):
        # Set SharPy output directory and file root.
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'dAgamma0_dzeta'


    def tearDown(self):
        pass
    
    def test_dAgamma0_dzeta_unit(self):
        """Test matrix against numerical example."""
        M=1
        N=1
        K=M*N
        dAgam0_dzeta = np.zeros((3*(M+1)*(N+1)))
        zeta=np.array((0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,1.0,1.0,0.0))
        gamma0=np.array((1.0))
        dZeta=np.random.random(len(zeta))/10.0 #10% variations
        # create matrix
        Cpp_dAgamma0_dZeta(zeta, M, N, gamma0, zeta, M, N, dAgam0_dzeta)
        # change in downwash
        dwApprox = np.dot(dAgam0_dzeta,dZeta)
        # calculate AICs for numrical solution
        AIC = np.zeros((K,K))
        AIC2= np.zeros((K,K))
        Cpp_AIC(zeta, M, N, zeta, M, N, AIC)
        Cpp_AIC(zeta+dZeta, M, N, zeta+dZeta, M, N, AIC2)
        # calculate exact solution
        dwExact = np.dot(AIC2,gamma0) - np.dot(AIC,gamma0)
        # check less than 1% maximum rel error
        self.assertLess(np.absolute((dwApprox - dwExact)/(np.dot(AIC,gamma0))),0.01)
        
    def test_dAgammaW0_dzeta_unit(self):
        """Test matrix against numerical example."""
        """Test matrix against numerical example."""
        M=1
        N=1
        K=M*N
        dAgamW0_dzeta = np.zeros((3*(M+1)*(N+1)))
        zeta=np.array((0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,1.0,1.0,0.0))
        zetaW=np.array((1.0,0.0,0.0,1.0,1.0,0.0,2.0,0.0,0.0,2.0,1.0,0.0))
        gammaW0=np.array((1.0))
        Cpp_dAgamma0_dZeta(zetaW, M, N, gammaW0, zeta, M, N, dAgamW0_dzeta)
        for foo in range(100):
            dZeta=np.random.random(len(zeta))/100.0 #1% variations
            # change in downwash
            dwApprox = np.dot(dAgamW0_dzeta,dZeta)
            # calculate AICs for numrical solution
            AIC = np.zeros((K,K))
            AIC2= np.zeros((K,K))
            Cpp_AIC(zetaW, M, N, zeta, M, N, AIC)
            Cpp_AIC(zetaW, M, N, zeta+dZeta, M, N, AIC2)
            # calculate exact solution
            dwExact = np.dot(AIC2,gammaW0) - np.dot(AIC,gammaW0)
            # check less than 1% maximum rel error
            if False:
                print("test no. ",foo)
                print(np.dot(AIC,gammaW0))
                print(dwExact)
                print(dwApprox)
                print((dwApprox - dwExact)/(np.dot(AIC,gammaW0)))
            # check less than 2%
            self.assertLess(np.absolute((dwApprox - dwExact)/(np.dot(AIC,gammaW0))),0.01)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()