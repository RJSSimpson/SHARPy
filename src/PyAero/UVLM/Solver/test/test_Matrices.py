'''
Created on 1 Jun 2015

@author: rjs10
'''
import unittest
import SharPySettings as Settings
import numpy as np
import ctypes as ct
from UVLMLib import Cpp_AIC
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

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()