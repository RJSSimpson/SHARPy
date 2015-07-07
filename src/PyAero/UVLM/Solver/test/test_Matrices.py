'''
Created on 1 Jun 2015

@author: rjs10
'''
import unittest
import SharPySettings as Settings
import numpy as np
import ctypes as ct
from UVLMLib import Cpp_AIC, Cpp_dAgamma0_dZeta, Cpp_genW, Cpp_dWzetaPri0_dZeta, Cpp_genH, Cpp_AIC3, Cpp_dA3gamma0_dZeta, Cpp_Y1
from scipy.io import savemat

TestDir = (Settings.SharPyProjectDir + 'output/tests/PyAero/UVLM/' 
           + 'Matrices/')

class Test_AIC(unittest.TestCase):
    """@brief Test AIC matrix elements.
       @notes The dAgamma0W_dZeta approximation seems to have more relative error
       than the dAgamma0_dZeta, but less than 1% for random 0-1% variations in every DoF."""

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
        
class Test_dAgamma0_dzeta(unittest.TestCase):
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
            
    def test_dAgamma0_dzeta_unit2DM10(self):
        """Test matrix against numerical example for 2D aerofoil, 10 panels."""
        #Init AIC
        m=10
        n=1
        k=m*n
        gamma0=np.ones((k))
        AIC = np.zeros((k,k))
        AIC2 = np.zeros((k,k))
        dAgam0_dzeta = np.zeros((k,3*(m+1)*(n+1)))
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
        # gen matrix
        Cpp_dAgamma0_dZeta(zeta, m, n, gamma0, zeta, m, n, dAgam0_dzeta)        
        for foo in range(10):
            # gen random zeta
            dZeta=np.random.random(len(zeta))/100.0 #10% chord panel vars
            # calc dw Approx
            dwApprox = np.dot(dAgam0_dzeta,dZeta)
            # calc AICs for numerical comparison
            Cpp_AIC(zeta, m, n, zeta, m, n, AIC)
            Cpp_AIC(zeta+dZeta, m, n, zeta+dZeta, m, n, AIC2)
            dwExact = np.dot(AIC2,gamma0) - np.dot(AIC,gamma0)
            if False:
                print("test no. ",foo)
                print(np.dot(AIC,gamma0))
                print(dwExact)
                print(dwApprox)
                print((dwExact-dwApprox)/np.dot(AIC,gamma0))
            self.assertLess(np.max(np.absolute((dwApprox - dwExact)/(np.dot(AIC,gamma0)))),0.01)
        
    def test_dAgammaW0_dzeta_unit2DM10(self):
        """Test matrix against numerical example for 2D aerofoil, 10 panels, 90 wake panels."""
        #Init AIC
        m=10
        mW=190
        n=2
        k=m*n
        kW=mW*n
        gammaW0=np.ones((kW))
        AIC = np.zeros((k,kW))
        AIC2 = np.zeros((k,kW))
        dAgamW0_dzeta = np.zeros((k,3*(m+1)*(n+1)))
        # initialize grid for 3D solver (Unit VR)
        chords = np.linspace(0.0, 1.0, m+1, True)
        chordsW = np.linspace(1.0,mW*(1.0/m),mW+1,True)
        spans = np.linspace(-1000, 1000, n+1,True)
        zeta=np.zeros(3*len(chords)*len(spans))
        zetaW=np.zeros(3*len(chordsW)*len(spans))
        kk=0
        for c in chords:
            for s in spans:
                zeta[3*kk]=c
                zeta[3*kk+1]=s
                kk=kk+1
        kk=0
        for c in chordsW:
            for s in spans:
                zetaW[3*kk]=c
                zetaW[3*kk+1]=s
                kk=kk+1
        # gen matrix
        Cpp_dAgamma0_dZeta(zetaW, mW, n, gammaW0, zeta, m, n, dAgamW0_dzeta)      
        for foo in range(10):
            # gen random zeta
            dZeta=np.random.random(len(zeta))/1000.0 #1% chord panel vars
            # calc dw Approx
            dwApprox = np.dot(dAgamW0_dzeta,dZeta)
            # calc AICs for numerical comparison
            Cpp_AIC(zetaW, mW, n, zeta, m, n, AIC)
            Cpp_AIC(zetaW, mW, n, zeta+dZeta, m, n, AIC2)
            dwExact = np.dot(AIC2,gammaW0) - np.dot(AIC,gammaW0)
            if False:
                print("test no. ",foo)
                print(np.dot(AIC,gammaW0))
                print(dwExact)
                print(dwApprox)
                print((dwExact-dwApprox)/np.dot(AIC,gammaW0))
            self.assertLess(np.max(np.absolute((dwApprox - dwExact)/(np.dot(AIC,gammaW0)))),0.01)


class Test_dWzetaPri0_dzeta(unittest.TestCase):
    """@brief Test non-zero reference zetaPri matrix, dWzetaPri0_dzeta."""
    
    def setUp(self):
        # Set SharPy output directory and file root.
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'dWzetaPri0_dzeta'

    def tearDown(self):
        pass
    
    def test_dWzetaPri0_dzeta_unit2DM10(self):
        """Test matrix for example of aerofoil, 10 chordwise panels."""
        m=10
        n=2
        k=m*n
        zetaPri0=np.ones((3*(m+1)*(n+1)))
        W = np.zeros((k,3*(m+1)*(n+1)))
        W[:,:]=0.0 # numpy not initializing to zero!!
        W2 = np.zeros((k,3*(m+1)*(n+1)))
        W2[:,:] = 0.0
        dWzetaPri0_dzeta = np.zeros((k,3*(m+1)*(n+1)))
        dWzetaPri0_dzeta[:,:] = 0.0
        # initialize grid for 3D solver (Unit VR)
        chords = np.linspace(0.0, 1.0, m+1, True)
        spans = np.linspace(-1000, 1000, n+1,True)
        zeta=np.zeros(3*len(chords)*len(spans))
        kk=0
        for c in chords:
            for s in spans:
                zeta[3*kk]=c
                zeta[3*kk+1]=s
                kk=kk+1
        # gen matrix
        Cpp_genW(zeta,m,n,W)
        Cpp_dWzetaPri0_dZeta(zeta, m, n, zetaPri0, dWzetaPri0_dzeta)
        for foo in range(10):
            # gen random zeta
            dZeta=np.random.random(len(zeta))/1000.0 #1% chord panel vars
            # calculate approx
            dwApprox = np.dot(dWzetaPri0_dzeta,dZeta)
            # claculate exact
            Cpp_genW(zeta+dZeta,m,n,W2)
            dwExact = np.dot(W2,zetaPri0) - np.dot(W,zetaPri0)
            if False:
                print("test no. ",foo)
                print(np.dot(W,zetaPri0))
                print(dwExact)
                print(dwApprox)
                print((dwExact-dwApprox)/np.dot(W,zetaPri0))
            self.assertLess(np.max(np.absolute((dwExact-dwApprox)/np.dot(W,zetaPri0))),0.01)
            
class Test_H(unittest.TestCase):
    """@brief Test segment to lattice vertex distribution matrix."""
    
    def setUp(self):
        # Set SharPy output directory and file root.
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'H'

    def tearDown(self):
        pass
    
    def test_HvKnown(self):
        """Test matrix against trusted result."""
        m=2
        n=1
        k=m*n
        kZeta=(m+1)*(n+1)
        H = np.zeros((3*kZeta,12*k))
        Cpp_genH(m,n,H)
        data = np.loadtxt(TestDir + "H_M2N1.dat")
        self.assertTrue(np.array_equal(H, data))
        
class Test_AIC3(unittest.TestCase):
    """@brief Test 3 component AIC."""
    
    def setUp(self):
        # Set SharPy output directory and file root.
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'AIC3'

    def tearDown(self):
        pass
    
    def test_AIC3vAIC(self):
        """Test against trusted result of typical AIC."""
        #Init AIC
        m=5
        n=1
        k=m*n
        AIC = np.zeros((k,k))
        AIC3 = np.zeros((3*k,k))
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
        Cpp_AIC3(zeta, m, n, zeta, m, n, AIC3)
        for i in range(k-1):
            for j in range(k-1):
                self.assertAlmostEqual(AIC3[3*i+2,j],AIC[i,j],4)
            # end for j
        # end for i
        
class Test_dA3gamma0_dZeta(unittest.TestCase):
    """@brief Test variations of 3 component AIC matrix."""
    
    def setUp(self):
        # Set SharPy output directory and file root.
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'dA3gamma0_dZeta'

    def tearDown(self):
        pass
    
    def test_dA3gamma0_dzeta_unit(self):
        """Test matrix against numerical example."""
        M=1
        N=1
        K=M*N
        dA3gam0_dzeta = np.zeros((3,3*(M+1)*(N+1)))
        zeta=np.array((0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,1.0,1.0,0.0))
        gamma0=np.array((1.0))
        dZeta=np.random.random(len(zeta))/10.0 #10% variations
        # create matrix
        Cpp_dA3gamma0_dZeta(zeta, M, N, gamma0, zeta, M, N, dA3gam0_dzeta)
        # change in downwash
        dwApprox = np.dot(dA3gam0_dzeta,dZeta)
        # calculate AICs for numrical solution
        AIC = np.zeros((3*K,K))
        AIC2= np.zeros((3*K,K))
        Cpp_AIC3(zeta, M, N, zeta, M, N, AIC)
        Cpp_AIC3(zeta+dZeta, M, N, zeta+dZeta, M, N, AIC2)
        # calculate exact solution
        dwExact = np.squeeze(np.dot(AIC2,gamma0) - np.dot(AIC,gamma0))
        # check less than 1% maximum rel error
        self.assertLess(np.max(np.absolute((dwApprox - dwExact)/np.linalg.norm(np.dot(AIC,gamma0)))),0.01)
        
    def test_dAgammaW0_dzeta_unit(self):
        """Test matrix against numerical example."""
        M=1
        N=1
        K=M*N
        dA3gamW0_dzeta = np.zeros((3,3*(M+1)*(N+1)))
        zeta=np.array((0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,1.0,1.0,0.0))
        zetaW=np.array((1.0,0.0,0.0,1.0,1.0,0.0,2.0,0.0,0.0,2.0,1.0,0.0))
        gammaW0=np.array((1.0))
        Cpp_dA3gamma0_dZeta(zetaW, M, N, gammaW0, zeta, M, N, dA3gamW0_dzeta)
        for foo in range(100):
            dZeta=np.random.random(len(zeta))/100.0 #1% variations
            # change in downwash
            dwApprox = np.dot(dA3gamW0_dzeta,dZeta)
            # calculate AICs for numrical solution
            AIC = np.zeros((3*K,K))
            AIC2= np.zeros((3*K,K))
            Cpp_AIC3(zetaW, M, N, zeta, M, N, AIC)
            Cpp_AIC3(zetaW, M, N, zeta+dZeta, M, N, AIC2)
            # calculate exact solution
            dwExact = np.squeeze(np.dot(AIC2,gammaW0) - np.dot(AIC,gammaW0))
            # check less than 1% maximum rel error
            if False:
                print("test no. ",foo)
                print(np.dot(AIC,gammaW0))
                print(dwExact)
                print(dwApprox)
                print((dwApprox - dwExact)/(np.linalg.norm(np.dot(AIC,gammaW0))))
            # check less than 2%
            self.assertLess(np.max(np.absolute((dwApprox - dwExact)/np.linalg.norm(np.dot(AIC,gammaW0)))),0.01)
            
    def test_dA3gamma0_dzeta_unit2DM10(self):
        """Test matrix against numerical example for 2D aerofoil, 10 panels."""
        #Init AIC
        m=10
        n=1
        k=m*n
        gamma0=np.ones((k))
        AIC = np.zeros((3*k,k))
        AIC2 = np.zeros((3*k,k))
        dA3gam0_dzeta = np.zeros((3*k,3*(m+1)*(n+1)))
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
        # gen matrix
        Cpp_dA3gamma0_dZeta(zeta, m, n, gamma0, zeta, m, n, dA3gam0_dzeta)        
        for foo in range(10):
            # gen random zeta
            dZeta=np.random.random(len(zeta))/100.0 #10% chord panel vars
            # calc dw Approx
            dwApprox = np.dot(dA3gam0_dzeta,dZeta)
            # calc AICs for numerical comparison
            Cpp_AIC3(zeta, m, n, zeta, m, n, AIC)
            Cpp_AIC3(zeta+dZeta, m, n, zeta+dZeta, m, n, AIC2)
            dwExact = np.dot(AIC2,gamma0) - np.dot(AIC,gamma0)
            if False:
                print("test no. ",foo)
                print(np.dot(AIC,gamma0))
                print(dwExact)
                print(dwApprox)
                print(dwExact-dwApprox)
            nuC = np.dot(AIC,gamma0)
            for kk in range(k):
                self.assertLess(np.max(np.absolute((dwExact-dwApprox)[3*kk:3*kk+3]/np.linalg.norm(nuC[3*kk:3*kk+3]))),0.01)

    def test_dA3gammaW0_dzeta_unit2DM10(self):
        """Test matrix against numerical example for 2D aerofoil, 10 panels, 90 wake panels."""
        #Init AIC
        m=10
        mW=190
        n=2
        k=m*n
        kW=mW*n
        gammaW0=np.ones((kW))
        AIC = np.zeros((3*k,kW))
        AIC2 = np.zeros((3*k,kW))
        dA3gamW0_dzeta = np.zeros((3*k,3*(m+1)*(n+1)))
        # initialize grid for 3D solver (Unit VR)
        chords = np.linspace(0.0, 1.0, m+1, True)
        chordsW = np.linspace(1.0,mW*(1.0/m),mW+1,True)
        spans = np.linspace(-1000, 1000, n+1,True)
        zeta=np.zeros(3*len(chords)*len(spans))
        zetaW=np.zeros(3*len(chordsW)*len(spans))
        kk=0
        for c in chords:
            for s in spans:
                zeta[3*kk]=c
                zeta[3*kk+1]=s
                kk=kk+1
        kk=0
        for c in chordsW:
            for s in spans:
                zetaW[3*kk]=c
                zetaW[3*kk+1]=s
                kk=kk+1
        # gen matrix
        Cpp_dA3gamma0_dZeta(zetaW, mW, n, gammaW0, zeta, m, n, dA3gamW0_dzeta)      
        for foo in range(10):
            # gen random zeta
            dZeta=np.random.random(len(zeta))/1000.0 #1% chord panel vars
            # calc dw Approx
            dwApprox = np.dot(dA3gamW0_dzeta,dZeta)
            # calc AICs for numerical comparison
            Cpp_AIC3(zetaW, mW, n, zeta, m, n, AIC)
            Cpp_AIC3(zetaW, mW, n, zeta+dZeta, m, n, AIC2)
            dwExact = np.dot(AIC2,gammaW0) - np.dot(AIC,gammaW0)
            if False:
                print("test no. ",foo)
                print(np.dot(AIC,gammaW0))
                print(dwExact)
                print(dwApprox)
            nuC = np.dot(AIC,gammaW0)
            for kk in range(k):
                self.assertLess(np.max(np.absolute((dwExact-dwApprox)[3*kk:3*kk+3]/np.linalg.norm(nuC[3*kk:3*kk+3]))),0.01)
                
class Test_Ys(unittest.TestCase):
    """@brief Test variations of 3 component AIC matrix."""
    
    def setUp(self):
        # Set SharPy output directory and file root.
        Settings.OutputDir = TestDir
        Settings.OutputFileRoot = 'Y'

    def tearDown(self):
        pass
    
    def test_Y1(self):
        """Test Y1 matrix generation."""
        # colloc velocities
        vC = np.array((0.0,0.0,1.0, 0.0,0.0,1.0, 0.0,0.0,1.0, 0.0,0.0,1.0))
        # zeta
        m=1
        n=1
        chords = np.linspace(0.0, 1.0, m+1, True)
        spans = np.linspace(0.0, 1.0, n+1, True)
        zeta=np.zeros(3*len(chords)*len(spans))
        kk=0
        for c in chords:
            for s in spans:
                zeta[3*kk]=c
                zeta[3*kk+1]=s
                kk=kk+1
        # init Y1
        Y1 = np.zeros((12*m*n,m*n))
        Cpp_Y1(vC, zeta, m, n, Y1)
        np.array_equal(Y1, [-1.0, 0., 0., 0., 1.0, 0.0, 1.0, 0., 0., 0., 1.0, 0.])

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    #unittest.main()
    suite1 = unittest.TestLoader().loadTestsFromTestCase(Test_Ys)
    alltests = unittest.TestSuite([suite1])
    TestRunner = unittest.TextTestRunner(verbosity=2)
    TestRunner.run(alltests)