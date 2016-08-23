'''
@package    PyAero.UVLM.UVLM
@brief      UVLM solution.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       18/02/2013
@pre        None
@warning    None
'''

import UVLMLib
import numpy as np
from DerivedTypesAero import VMopts, VMinput, ControlSurf
import ctypes as ct
from XbeamLib import Psi2TransMat
import DerivedTypes
import PostProcess
import SharPySettings as Settings
from VLM import InitSteadyExternalVels, InitSteadyWake
import BeamInit
from PyFSI.Beam2UVLM import InitSection, CoincidentGrid
from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps
from PyAero.UVLM.Utils.DerivedTypesAero import VMUnsteadyInput
from scipy.io import savemat, loadmat
import getpass
from PyBeam.Utils.Misc import iNode2iElem

def Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS,vOmegaHist=None,eta0=None,numNodesElem=2,delEtaHist=None):
    """@brief UVLM solver with prescribed inputs.
       @param vOmegaHist History of velocities of the A-frame, expressed in A-frame.
       @param eta0 History of underlying beam disps/rotations in A-frame."""
    
    # Initialise DerivedTypes for PyBeam initialisation. 
    XBOPTS = DerivedTypes.Xbopts()
    if numNodesElem==2:
        beamElems=VMOPTS.N.value
    elif numNodesElem==3:
        beamElems=int(VMOPTS.N.value/2)
    else:
        raise ValueError('numNodesElem should be 2 or 3')
        
    XBINPUT = DerivedTypes.Xbinput(numNodesElem, beamElems, VMINPUT.b)  
    # Create a dummy stiffness matrix to avoid errors.
    XBINPUT.BeamStiffness = np.eye(6,6,0,ct.c_double)
    # Use PyBeam to define reference beam (elastic axis).
    XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni,\
            XBNODE, NumDof = BeamInit.Static(XBINPUT, XBOPTS)
    # Delete unnecessary variables.
    del XBELEM, XBNODE, NumDof  
    # Copy reference beam to current (deformed) variables.
    if eta0==None:
        PosDefor = PosIni.copy(order = 'F')
        PsiDefor = PsiIni.copy(order = 'F')
        # Declare empty array for beam DoF rates.
        PosDotDef = np.zeros_like(PosDefor, ct.c_double, 'F')
        PsiDotDef = np.zeros_like(PsiDefor, ct.c_double, 'F')
    else:
        PosDefor=eta0[0]+delEtaHist[0][:,:,0]
        PsiDefor=eta0[1]+delEtaHist[1][:,:,:,0]
        PosDotDef=delEtaHist[2][:,:,0]
        PsiDotDef=delEtaHist[3][:,:,:,0]
       
    # Check the specified inputs from the PyAero have been properly applied.
    assert NumNodes_tot.value - 1 == VMOPTS.N.value, "Initialisation wrong"
    
    # Initialise section coordinates.
    Section = InitSection(VMOPTS,VMINPUT,AELOPTS.ElasticAxis)
    
    # Initialise origin and orientation of surface velocities in a-frame
    CGa = Psi2TransMat(VMUNST.PsiA_G)
    if vOmegaHist == None:
        VelA_A = np.dot(CGa.T,VMUNST.VelA_G)
        OmegaA_A = np.dot(CGa.T,VMUNST.OmegaA_G)
    else:
        VelA_A = np.dot(CGa.T,vOmegaHist[0,1:4])
        OmegaA_A = np.dot(CGa.T,vOmegaHist[0,4:7])
    
    # Declare empty array for aerodynamic grid and velocities.
    Zeta = np.zeros((VMOPTS.M.value + 1,VMOPTS.N.value + 1,3), ct.c_double, 'C')
    ZetaDot = np.zeros((VMOPTS.M.value + 1, VMOPTS.N.value + 1, 3),
                       ct.c_double, 'C')
    
    # Initialise aerodynamic grid and velocities.
    CoincidentGrid(PosDefor, PsiDefor, Section,
                   VelA_A, OmegaA_A, PosDotDef, PsiDotDef,
                   XBINPUT, Zeta, ZetaDot, VMUNST.OriginA_G, VMUNST.PsiA_G,
                   VMINPUT.ctrlSurf)
    
    # Initialise wake for unsteady solution.
    if vOmegaHist == None:
        velForWake = VMUNST.VelA_G
    else:
        velForWake = vOmegaHist[0,1:4]
        
    ZetaStar, GammaStar = InitSteadyWake(VMOPTS, VMINPUT, Zeta, velForWake)
    
    # Initialise external velocities.  
    Uext = InitSteadyExternalVels(VMOPTS,VMINPUT)
    
    # Declare empty solver variables.
    Forces = np.zeros_like(Zeta, ct.c_double, 'C')
    Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')
    AIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,
                   VMOPTS.M.value*VMOPTS.N.value),
                   ct.c_double,'C')
    BIC = np.zeros((VMOPTS.M.value*VMOPTS.N.value,
                   VMOPTS.M.value*VMOPTS.N.value),
                   ct.c_double,'C')
    
    # Open tecplot file object."
    Variables=['X', 'Y', 'Z','Gamma']
    Filename = Settings.OutputDir + Settings.OutputFileRoot + 'AeroGrid.dat'
    FileObject = PostProcess.WriteAeroTecHeader(Filename,
                                                'Default',
                                                Variables)
    
    # Initialise vector of time steps.
    Time = np.arange(0.0,VMUNST.FinalTime,VMOPTS.DelTime.value)
    # Create Array for storing time and coefficient data.
    CoeffHistory = np.zeros((len(Time),4))
    
    # Loop through time steps.
    for iTimeStep in range(len(Time)):
        
        # Set forces array to zero (+= operator used in C++ library).
        Forces[:,:,:] = 0.0
        
        if iTimeStep > 0:
            
            # Update geometry.
            if vOmegaHist == None:
                VMUNST.OriginA_G[:] += VMUNST.VelA_G[:]*VMUNST.DelTime
                # TODO: update OmegaA_A, PsiA_A in pitching problem.
            else:
                VMUNST.OriginA_G[:] += vOmegaHist[iTimeStep,1:4]*VMUNST.DelTime
                VelA_A = vOmegaHist[iTimeStep,1:4]
                # TODO: update OmegaA_A, PsiA_G in pitching problem.
                
            if delEtaHist != None:
                PosDefor=eta0[0]+delEtaHist[0][:,:,iTimeStep]
                PsiDefor=eta0[1]+delEtaHist[1][:,:,:,iTimeStep]
                PosDotDef=delEtaHist[2][:,:,iTimeStep]
                PsiDotDef=delEtaHist[3][:,:,:,iTimeStep]
            
            # Update control surface defintion.
            if VMINPUT.ctrlSurf != None:
                VMINPUT.ctrlSurf.update(Time[iTimeStep])
            
            # Update aerodynamic surface.
            CoincidentGrid(PosDefor,
                           PsiDefor,
                           Section,
                           VelA_A,
                           OmegaA_A,
                           PosDotDef,
                           PsiDotDef,
                           XBINPUT,
                           Zeta,
                           ZetaDot,
                           VMUNST.OriginA_G,
                           VMUNST.PsiA_G,
                           VMINPUT.ctrlSurf)
            
            # Convect wake downstream.        
            ZetaStar = np.roll(ZetaStar,1,axis = 0)
            GammaStar = np.roll(GammaStar,1,axis = 0)
            # Overwrite 1st row with with new trailing-edge position.
            ZetaStar[0,:] = Zeta[VMOPTS.M.value,:]
            # Overwrite Gamma with TE value from previous time step.
            GammaStar[0,:] = Gamma[VMOPTS.M.value-1,:]
            
            # set NewAIC to false.
            if VMINPUT.ctrlSurf == False:
                VMOPTS.NewAIC = ct.c_bool(False)
            
        # END if iTimeStep > 1
        
        
        # Calculate forces on aerodynamic grid.
        if iTimeStep == 0:
            sav = VMOPTS.NewAIC
            VMOPTS.NewAIC = ct.c_bool(1)
        
        UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS,
                               Forces, Gamma, GammaStar, AIC, BIC)
        
        if iTimeStep == 0:
            VMOPTS.NewAIC = sav
            del sav
         
        #print(PostProcess.GetCoeffs(VMOPTS, Forces, VMINPUT, VMUNST.VelA_G))
        
        CoeffHistory[iTimeStep,0] = Time[iTimeStep]
        CoeffHistory[iTimeStep,1:] = PostProcess.GetCoeffs(VMOPTS, Forces,
                                                           VMINPUT,
                                                           VMUNST.VelA_G)
        
        # Write aerodynamic surface data as tecplot zone data.
        PostProcess.WriteAeroTecZone(FileObject, 'Surface', Zeta,\
                                     iTimeStep, len(Time),\
                                     Time[iTimeStep], Variables, False, Gamma)
        
        # Write wake data as tecplot zone data.
        PostProcess.WriteAeroTecZone(FileObject, 'Wake', ZetaStar,\
                                    iTimeStep, len(Time),\
                                    Time[iTimeStep], \
                                    Variables, False, GammaStar)
        
        # Rollup due to external velocities.
        ZetaStar[:,:] += VMINPUT.U_infty * VMOPTS.DelTime.value
        
    # END for iTimeStep
    
    # Close tecplot file object.
    PostProcess.CloseAeroTecFile(FileObject)
    
    # write geometry and circ dist. to file for matlab
    if iTimeStep == len(Time)-1 and False:
        savemat('/home/rjs10/Desktop/uvlmOut.mat',
                {'zeta':Zeta.flatten('C'),
                 'zetaW':ZetaStar.flatten('C'),
                 'gamma':Gamma.flatten('C'),
                 'gammaW':GammaStar.flatten('C')},
                False,
                oned_as='column'
               )
    
    return CoeffHistory
        

if __name__ == '__main__':
    Settings.OutputDir = '/home/' + getpass.getuser() + '/Documents/MATLAB/Patil_HALE/nonZeroAero/'
    writeToMat = True
    
    # Define options for aero solver.
    M = 4
    N = 10
    Mstar = 10*M
    U_mag = 25.0
    alpha = 2.0*np.pi/180.0
    theta = 0.0*np.pi/180.0
    imageMeth=True
    
    # Define wing geometry parameters.
    c = 1
    b = 16 #semi-span
    
    # Physical time-step
    DeltaS = c/(M*U_mag)
    
    # Solver options
    VMOPTS = VMopts(M,N,imageMeth,Mstar,False,True,False,DeltaS,False,4)
    
    # Solver inputs
    WakeLength = Mstar/M # specified in chord lengths
    ctrlSurf = None #ControlSurf(3, 4, 6, 9, 'sin', 1*np.pi/180.0, 2*np.pi)
    VMINPUT = VMinput(c, b, U_mag, alpha, theta, WakeLength,  ctrlSurf)
    
    
    # Define unsteady solver parameters.
    NumChordLengths = 200.0
    VelA_A = np.array([0, 0, 0])
    OmegaA_A = np.array([0, 0, 0])
    VMUNST = VMUnsteadyInput(VMOPTS,VMINPUT,\
                             WakeLength,\
                             DeltaS*U_mag/c,\
                             NumChordLengths,\
                             VelA_A, OmegaA_A)
    
    # Generate history of axis sytem motions (in inertial frame)
    Time = np.arange(0.0,VMUNST.FinalTime,VMOPTS.DelTime.value)
    if False:
        hBar=0.01
        k=1.0
        omegaY=2*U_mag*k/c
        vOmegaHist = np.zeros((len(Time),7))
        vOmegaHist[:,0] = Time
        vOmegaHist[:,1] = 0.0 # along A-frame x-axis
        vOmegaHist[:,2] = omegaY*hBar*np.cos(omegaY*Time) # y-axis (fore/aft)
        vOmegaHist[:,3] = 0.0 # z-axis (up/down)
        vOmegaHist[:,4] = 0.0 # about A-frame x-axis
        vOmegaHist[:,5] = 0.0 # about y
        vOmegaHist[:,6] = 0.0 # about z
    else:
        vOmegaHist = np.zeros((len(Time),7))
        vOmegaHist[:,0] = Time
        
    if True:
        refFile = '/home/rjs10/Documents/MATLAB/Patil_HALE/nonZeroAero/M4N10_V25_alpha2_SOL112_def.dat'
        fp = open(refFile)
        numNodesElem = 3
        if numNodesElem == 2:
            NumElems = N
        elif numNodesElem == 3:
            NumElems = int(N/2)
        else:
            raise ValueError('numNodesElem should be 2 or 3')
        
        # Initialize PosDefor, PsiDefor
        PosDefor = np.zeros((N+1,3),dtype=ct.c_double, order='F')
        PosDeforElem = np.zeros((NumElems,3,3),dtype=ct.c_double, order='F')
        PsiDefor = np.zeros((NumElems,3,3),dtype=ct.c_double, order='F')
        
        data = np.loadtxt(refFile, skiprows=2)
        i0=0
        for iElem in range(NumElems):
            for i in range(numNodesElem):
                PosDeforElem[iElem,i,:] = data[i0,2:5]
                PsiDefor[iElem,i,:] = data[i0,5:8]
                i0+=1
        
        # extract nodal information (PosDefor)
        for iNode in range(N+1):
            iElem, iiElem = iNode2iElem(iNode, N+1, numNodesElem)
            PosDefor[iNode,:] = PosDeforElem[iElem,iiElem,:]
        
        eta0 = (PosDefor, PsiDefor)
    else:
        eta0 = None
        
    # Generate a history of beam DoF motions (in A-frame)
    if True:
        delEtaHisty = np.zeros((12*N,len(Time)))
        PosDefHist=np.zeros((N+1,3,len(Time)),dtype=ct.c_double, order='F')
        PosDotHist=np.zeros((N+1,3,len(Time)),dtype=ct.c_double, order='F')
        PsiDefHistNode=np.zeros((N+1,3,len(Time)),dtype=ct.c_double, order='F')
        PsiDotHistNode=np.zeros((N+1,3,len(Time)),dtype=ct.c_double, order='F')
        PsiDeforHist = np.zeros((NumElems,3,3,len(Time)),dtype=ct.c_double, order='F')
        PsiDotHist = np.zeros((NumElems,3,3,len(Time)),dtype=ct.c_double, order='F')
        # load and scale eigenvector
        modeFile = '/home/rjs10/Documents/MATLAB/Patil_HALE/nonZeroAero/M4N10_V25_alpha2_SSbeam'
        beamDict = loadmat(modeFile)
        modeNum = 1
        eigVec = beamDict['phiSort'][:,modeNum-1]
        disps=np.array((np.arange(0, np.size(eigVec, 0)-1, 6),
                        np.arange(1, np.size(eigVec, 0)-1, 6),
                        np.arange(2, np.size(eigVec, 0)-1, 6)))
        disps = np.sort(disps.flatten())
        rots = np.setdiff1d(range(6*N), disps)
        maxDef = np.max(eigVec[disps])
        scale=(1/maxDef)*0.01*16.0
        phi0=scale*eigVec
        # get mode freq
        k = np.real(np.sqrt(beamDict['kMat'][modeNum-1,modeNum-1]))
        omega = 2*U_mag*k/c
        # create history
        tCount=0
        for t in Time:
            delEtaHisty[:6*N,tCount]=phi0*np.sin(omega*t)
            delEtaHisty[6*N:,tCount]=omega*phi0*np.cos(omega*t)
            tCount+=1
        
        tCount=0
        for t in Time:
            for j in range(N):
                # displacements are saved nodally
                PosDefHist[j+1,:,tCount] = delEtaHisty[disps[3*j:3*j+3],tCount]
                PosDotHist[j+1,:,tCount] = delEtaHisty[6*N+disps[3*j:3*j+3],tCount]
                PsiDefHistNode[j+1,:,tCount] = delEtaHisty[rots[3*j:3*j+3],tCount]
                PsiDotHistNode[j+1,:,tCount] = delEtaHisty[6*N+rots[3*j:3*j+3],tCount]
            # end for j
            tCount+=1
        # end for t
        
        # populate connectivity array
        conn = np.zeros((Settings.MaxElNod*NumElems))
        if numNodesElem == 2:
            for ElemNo in range(0,NumElems):
                i = ElemNo*Settings.MaxElNod
                conn[i]=ElemNo+1
                conn[i+1]=ElemNo+2
                
            NumNodes_tot = ct.c_int(NumElems + 1) 
                
        elif numNodesElem == 3:
            for ElemNo in range(0,NumElems):
                i = ElemNo*Settings.MaxElNod
                conn[i]=2*ElemNo+1
                conn[i+1]=2*ElemNo+3
                conn[i+2]=2*ElemNo+2
                
        # loop through XBELEM.Conn to assign to [iElem,iiElem,:,tCount]
        for tCount in range(len(Time)):
            connCount=0
            for iElem in range(NumElems):
                for iiElem in range(numNodesElem):
                    PsiDeforHist[iElem,iiElem,:,tCount]=PsiDefHistNode[int(conn[connCount]-1),:,tCount]
                    PsiDotHist[iElem,iiElem,:,tCount]=PsiDotHistNode[int(conn[connCount]-1),:,tCount]
                    connCount+=1
        
        delEtaHist = (PosDefHist,PsiDeforHist,PosDotHist,PsiDotHist)
    
    else:
        delEtaHist = None
        
    # Define 'aeroelastic' options.
    AELOPTS = AeroelasticOps(AirDensity = 0.08891)
    
    # Run C++ solver
    Coeffs = Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT,VMUNST,AELOPTS,vOmegaHist,eta0,numNodesElem,delEtaHist)
    
    print(Coeffs[-50:,:])
    
    if writeToMat == True:
        fileName = Settings.OutputDir + 'UVLMrectAR' + str(b/c) + '_m' + str(M) + 'mW' + str(Mstar) + 'n' + str(N) + 'delS' + str(2/M) + 'V' + str(U_mag) + '_alpha' + str(alpha)[:6]
        if imageMeth != False:
            fileName += '_half'
        savemat(fileName,
                {'vOmegaHist':vOmegaHist,
                 'Coeffs':Coeffs},
                True)
    