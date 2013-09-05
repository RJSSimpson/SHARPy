'''@package PyCoupled.Coupled_NlnStatic
@brief      NonlinearStatic Beam + UVLM.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       11/02/2013
@pre        None
@warning    None
'''

import sys
import SharPySettings as Settings
import DerivedTypes
import BeamIO
import BeamLib
import BeamInit
import numpy as np
import ctypes as ct
from PyFSI.Beam2UVLM import InitSection
from PyFSI.Beam2UVLM import CoincidentGrid
from PyAero.UVLM.Utils import UVLMLib
from PyFSI.Beam2UVLM import CoincidentGridForce
from PyAero.UVLM.Utils import DerivedTypesAero
from PyAero.UVLM.Utils import PostProcess
from PyAero.UVLM.Solver.VLM import InitSteadyExternalVels
from PyAero.UVLM.Solver.VLM import InitSteadyWake
from PyCoupled.Utils.DerivedTypesAeroelastic import AeroelasticOps
from PyBeam.Utils.XbeamLib import Psi2TransMat

def isodd(num):
    """@brief returns True if Num is odd"""
    return bool(num & 1)

#TODO: Move to Beam Utils.
def AddGravityLoads(BeamForces,XBINPUT,XBELEM,AELAOPTS,PsiDefor,
                      chord = 1.0):
    """@brief Apply nodal gravity loads.
    @param BeamForces Nodal forces to update.
    @param XBINPUT Xbeam inputs.
    @param XBELEM Xbeam element information.
    @param AELAOPTS Aeroelastic options.
    @param PsiA_G Cartesian rotation vector describing orientation of a-frame
           with respect to earth.
    @param chord Aerofoil chord, assumed to be 1.0 if ommited. 
    @details Offset of mass centroid from elastic axis is currently calculated
             using the AELAOPTS.ElasticAxis and .InertialAxis parameters which 
             are analogous to that used by Theodorsen. It therefore assumes the 
             aerofoil section is defined on the y-axis and hence the moment arm
             points in the B-frame y-direction.
    @warning Assumes even distribution of nodes along beam.
    @warning Not valid for twisted aerofoil sections, i.e those that are defined
             with some theta angle.
    """
    
    MassPerLength = XBINPUT.BeamMass[0,0]
    
    if XBINPUT.NumNodesElem == 2:
        NodeSpacing = XBELEM.Length[0]
    elif XBINPUT.NumNodesElem == 3:
        NodeSpacing = 0.5*XBELEM.Length[0]
        
    ForcePerNode = -NodeSpacing * MassPerLength * XBINPUT.g 
    
    # Obtain transformation from Earth to a-frame.
    CGa = Psi2TransMat(XBINPUT.PsiA_G)
    CaG = CGa.T
    
    # Force in a-frame.
    Force_a = np.dot(CaG,np.array([0.0,0.0,ForcePerNode]))
    
    # Indices for boundary nodes.
    Root = 0
    Tip = BeamForces.shape[0]-1
    
    # Apply forces.
    BeamForces[Root+1:Tip,:3] += Force_a
    BeamForces[Root,:3] += 0.5*Force_a
    BeamForces[Tip,:3] += 0.5*Force_a
    
    # Get number of nodes per beam element.
    NumNodesElem = XBINPUT.NumNodesElem
    
    # Loop through nodes to get moment arm at each.
    for iNode in range(XBINPUT.NumNodesTot):
        
        # Work out what element we are in (works for 2 and 3-noded).
        if iNode == 0:
            iElem = 0
        elif iNode < XBINPUT.NumNodesTot-1:
            iElem = int(iNode/(NumNodesElem-1))
        elif iNode == XBINPUT.NumNodesTot-1:
            iElem = int((iNode-1)/(NumNodesElem-1))
            
        # Work out what sub-element node (iiElem) we are in.
        if NumNodesElem == 2:
            if iNode < XBINPUT.NumNodesTot-1:
                iiElem = 0
            elif iNode == XBINPUT.NumNodesTot-1:
                iiElem = 1
                
        elif NumNodesElem == 3:
            iiElem = 0
            if iNode == XBINPUT.NumNodesTot-1:
                iiElem = 2 
            elif isodd(iNode):
                iiElem = 1
        
        # Calculate transformation matrix for each node.
        CaB = Psi2TransMat(PsiDefor[iElem,iiElem,:])
        
        # Define moment arm in B-frame
        # Moment arm to elastic axis defined using EA and IA of Theodorsen.
        armY = -(AELAOPTS.InertialAxis - AELAOPTS.ElasticAxis)*chord/2.0
        armY_a = np.dot(CaB,np.array([0.0, armY, 0.0]))
        
        # Calculate moment
        if (iNode == Root or iNode == Tip):
            BeamForces[iNode,3:] += np.cross(armY_a, 0.5*Force_a)
        else:
            BeamForces[iNode,3:] += np.cross(armY_a, Force_a)
    

def Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,AELAOPTS):
        """Nonlinear static solver using Python to solve aeroelastic
        equation. Assembly of structural matrices is carried out with 
        Fortran subroutines. Aerodynamics solved using PyAero.UVLM."""
        
        
        """BEAM INIT --------------------------------------------------------"""
        
        "Check correct solution code"
        assert XBOPTS.Solution.value == 112, ('NonlinearStatic requested' +\
                                                  ' with wrong solution code')

        
        "Initialise beam"
        XBINPUT, XBOPTS, NumNodes_tot, XBELEM, PosIni, PsiIni, XBNODE, NumDof \
                    = BeamInit.Static(XBINPUT,XBOPTS)
                    
        
        "Set initial conditions as undef config"
        PosDefor = PosIni.copy(order='F')
        PsiDefor = PsiIni.copy(order='F')
        
        
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('Solve nonlinear static case in Python ... \n')
        
        "Initialise structural eqn tensors"
        KglobalFull = np.zeros((NumDof.value,NumDof.value),\
                                ct.c_double, 'F'); ks = ct.c_int()
        FglobalFull = np.zeros((NumDof.value,NumDof.value),\
                                ct.c_double, 'F'); fs = ct.c_int()
        
        DeltaS  = np.zeros(NumDof.value, ct.c_double, 'F')
        Qglobal = np.zeros(NumDof.value, ct.c_double, 'F')
        x       = np.zeros(NumDof.value, ct.c_double, 'F')
        dxdt    = np.zeros(NumDof.value, ct.c_double, 'F')
        
        "Beam Load Step tensors"
        iForceStep     = np.zeros((NumNodes_tot.value,6), ct.c_double, 'F')
        iForceStep_Dof = np.zeros(NumDof.value, ct.c_double, 'F')
               
        """END BEAM INIT ----------------------------------------------------"""
        
        
        """AERO INIT --------------------------------------------------------"""
 
        "define section used to generate aero mesh"       
        Section = InitSection(VMOPTS,VMINPUT,AELAOPTS.ElasticAxis)
     
        "declare memory for Aero grid and velocities"
        Zeta = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')
        ZetaDot = np.zeros((Section.shape[0],PosDefor.shape[0],3),ct.c_double,'C')
        
        "Additional Aero solver variables"
        AeroForces = np.zeros((VMOPTS.M.value+1,VMOPTS.N.value+1,3),ct.c_double,'C')
        Gamma = np.zeros((VMOPTS.M.value,VMOPTS.N.value),ct.c_double,'C')

        "init external velocities"  
        Uext = InitSteadyExternalVels(VMOPTS,VMINPUT)
        
        "create zero triads for motion of reference frame"
        VelA_A = np.zeros((3))
        OmegaA_A = np.zeros((3))
        
        "create zero vectors for structural vars not used in static analysis"
        PosDotDef = np.zeros_like(PosDefor, ct.c_double, 'F')
        PsiDotDef = np.zeros_like(PsiDefor, ct.c_double, 'F')      
        
        "define tecplot stuff"
        FileName = Settings.OutputDir + Settings.OutputFileRoot + 'AeroGrid.dat'
        Variables = ['X', 'Y', 'Z','Gamma']        
        FileObject = PostProcess.WriteAeroTecHeader(FileName, \
                                                    'Default',\
                                                    Variables)
    
        
        "Start Load Loop"
        for iLoadStep in range(XBOPTS.NumLoadSteps.value):
            
            "Reset convergence parameters"
            Iter = 0
            ResLog10 = 0.0
            
            "zero loads"
            XBINPUT.ForceStatic[:,:] = 0.0
            AeroForces[:,:,:] = 0.0
            
            """"calculate aero loads"""
            
            "create new grid"
            CoincidentGrid(PosDefor, PsiDefor, Section, VelA_A, \
                       OmegaA_A, PosDotDef, PsiDotDef, XBINPUT,\
                       Zeta, ZetaDot)
          
            "init wake"
            ZetaStar, GammaStar = InitSteadyWake(VMOPTS,VMINPUT,Zeta)
                 
            "solve for AeroForces"
            UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, \
                           AeroForces, Gamma, GammaStar)
            
            "apply density scaling"
            AeroForces[:,:,:] = AELAOPTS.AirDensity*AeroForces[:,:,:]
            
            "write to tecplot file"
            "surface"
            PostProcess.WriteAeroTecZone(FileObject, 'Surface', Zeta,\
                              iLoadStep,\
                              XBOPTS.NumLoadSteps.value,\
                              iLoadStep*1.0,\
                              Variables, True, Gamma)
            "wake"
            PostProcess.WriteAeroTecZone(FileObject, 'Wake', ZetaStar,\
                              iLoadStep,\
                              XBOPTS.NumLoadSteps.value,\
                              iLoadStep*1.0,\
                              Variables, False, GammaStar)

            
            # Map AeroForces to beam.
            CoincidentGridForce(XBINPUT, PsiDefor, Section, AeroForces,
                                XBINPUT.ForceStatic)
            
            # Add gravity loads.
            AddGravityLoads(XBINPUT.ForceStatic,XBINPUT,XBELEM,AELAOPTS,
                            PsiDefor,VMINPUT.c)

            
            "apply factor corresponding to force step"
            iForceStep = XBINPUT.ForceStatic*float( (iLoadStep+1) ) / \
                                                XBOPTS.NumLoadSteps.value
            
            if XBOPTS.PrintInfo.value == True:
                sys.stdout.write('  iLoad: %-10d\n' %(iLoadStep+1))
                sys.stdout.write('   SubIter DeltaF     DeltaX     ResLog10\n')
                
            "Newton Iteration"
            while( (ResLog10 > np.log10(XBOPTS.MinDelta.value)) \
                 & (Iter < XBOPTS.MaxIterations.value) ):
                
                "Increment iteration counter"
                Iter += 1
                if XBOPTS.PrintInfo.value == True:
                    sys.stdout.write('   %-7d ' %(Iter))
                    
                
                "Set structural eqn tensors to zero"
                KglobalFull[:,:] = 0.0; ks = ct.c_int()
                FglobalFull[:,:] = 0.0; fs = ct.c_int()
                Qglobal[:] = 0.0
                
                
                "Assemble matrices for static problem"
                BeamLib.Cbeam3_Asbly_Static(XBINPUT, NumNodes_tot, XBELEM, XBNODE,\
                            PosIni, PsiIni, PosDefor, PsiDefor,\
                            iForceStep, NumDof,\
                            ks, KglobalFull, fs, FglobalFull, Qglobal,\
                            XBOPTS)
                
                
                "Get state vector from current deformation"
                PosDot = np.zeros((NumNodes_tot.value,3), ct.c_double, 'F') #R
                PsiDot = np.zeros((XBINPUT.NumElems,Settings.MaxElNod,3),\
                                   ct.c_double, 'F') #Psi
                
                BeamLib.Cbeam_Solv_Disp2State(NumNodes_tot, NumDof, XBINPUT, XBNODE,\
                              PosDefor, PsiDefor, PosDot, PsiDot,
                              x, dxdt)
                
                
                "Get forces on unconstrained nodes"
                BeamLib.f_fem_m2v(ct.byref(NumNodes_tot),\
                                  ct.byref(ct.c_int(6)),\
                                  iForceStep.ctypes.data_as(ct.POINTER(ct.c_double)),\
                                  ct.byref(NumDof),\
                                  iForceStep_Dof.ctypes.data_as(ct.POINTER(ct.c_double)),\
                                  XBNODE.Vdof.ctypes.data_as(ct.POINTER(ct.c_int)) )
                
                "Calculate \Delta RHS"
                Qglobal = Qglobal - np.dot(FglobalFull, iForceStep_Dof)
                
                
                "Calculate \Delta State Vector"    
                DeltaS = - np.dot(np.linalg.inv(KglobalFull), Qglobal)
                    
                if XBOPTS.PrintInfo.value == True:
                    sys.stdout.write('%-10.4e %-10.4e ' \
                                     % (max(abs(Qglobal)),max(abs(DeltaS))))
                
                
                "Update Solution"
                BeamLib.Cbeam3_Solv_Update_Static(XBINPUT, NumNodes_tot, XBELEM,\
                                                  XBNODE, NumDof, DeltaS,\
                                                  PosIni, PsiIni, PosDefor,PsiDefor)
                
                
                "Residual at first iteration"
                if(Iter == 1):
                    Res0_Qglobal = max(abs(Qglobal)) + 1.e-16
                    Res0_DeltaX  = max(abs(DeltaS)) + 1.e-16
    
                "Update residual and compute log10"
                Res_Qglobal = max(abs(Qglobal)) + 1.e-16
                Res_DeltaX  = max(abs(DeltaS)) + 1.e-16
    
                ResLog10 = max([np.log10(Res_Qglobal/Res0_Qglobal), \
                                np.log10(Res_DeltaX/Res0_DeltaX)])
                
                
                if XBOPTS.PrintInfo.value == True:
                    sys.stdout.write('%8.4f\n' %(ResLog10))
                    
                "Stop the solution"
                if(ResLog10 > 10.):
                    sys.stderr.write(' STOP\n')
                    sys.stderr.write(' The max residual is %e\n' %(ResLog10))
                    exit(1)
                
            "END Newton Loop"
        "END Load Loop"
    
    
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write(' ... done\n')
        
         
        "Write deformed configuration to file"
        ofile = Settings.OutputDir + Settings.OutputFileRoot + '_SOL112_def.dat'
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('Writing file %s ... ' %(ofile))
        fp = open(ofile,'w')
        fp.write('TITLE="Non-linear static solution: deformed geometry"\n')
        fp.write('VARIABLES="iElem" "iNode" "Px" "Py" "Pz" "Rx" "Ry" "Rz"\n')
        fp.close()
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('done\n')
        WriteMode = 'a'
        
        BeamIO.OutputElems(XBINPUT.NumElems, NumNodes_tot.value, XBELEM, \
                           PosDefor, PsiDefor, ofile, WriteMode)
        
        
        "Print deformed configuration"
        if XBOPTS.PrintInfo.value==True:
            sys.stdout.write('--------------------------------------\n')
            sys.stdout.write('NONLINEAR STATIC SOLUTION\n')
            sys.stdout.write('%10s %10s %10s\n' %('X','Y','Z'))
            for inodi in range(NumNodes_tot.value):
                sys.stdout.write(' ')
                for inodj in range(3):
                    sys.stdout.write('%12.5e' %(PosDefor[inodi,inodj]))
                sys.stdout.write('\n')
            sys.stdout.write('--------------------------------------\n')
            
        "Close Tecplot ascii FileObject"
        PostProcess.CloseAeroTecFile(FileObject)
            
        
        "Return solution as optional output argument"
        return PosDefor, PsiDefor, Zeta, ZetaStar, Gamma, GammaStar, iForceStep

if __name__ == '__main__':
    "solve nlnstatic problem"
    
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
                                     KJMeth = True)
    
    "aero inputs"
    VMINPUT = DerivedTypesAero.VMinput(c = 1.0, b = XBINPUT.BeamLength,\
                                       U_mag = 25.0,\
                                       alpha = 4.0*np.pi/180.0,\
                                       theta = 0.0)
    
    "aeroelastic opts"
    "Density and gravity due to US standard atmosphere at 20km"
    AELAOPTS = AeroelasticOps(0.0,0.0,0.08891)
    
    
    Solve_Py(XBINPUT,XBOPTS,VMOPTS,VMINPUT,AELAOPTS)
    
    
    