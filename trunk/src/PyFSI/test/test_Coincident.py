'''@package PyFSI.Beam2UVLM.test_Coincident
@brief      print out deformed aero grid.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       21/01/2013
@pre        None
@warning    None
'''
import Beam2UVLM
from SectionInit import FlatPlateReg
import DerivedTypes
import NonlinearStatic
import numpy as np
import SharPySettings as Settings

if __name__ == '__main__':
    """Set up Xbopts for nonlinear static analysis defined in 
    NonlinearStatic/testcases.pdf case 1.1"""
    XBOPTS = DerivedTypes.Xbopts()
    XBOPTS.Solution.value = 112 
    XBOPTS.NumLoadSteps.value = 10
    XBOPTS.MinDelta.value = 1e-05
    XBOPTS.FollowerForce.value = False
         
    """Set up Xbinput for nonlinear static analysis defined in 
    NonlinearStatic/testcases.pdf case 1.1 with extra out-of plane loads."""
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
    XBINPUT.BeamMass[4,4] = 0.0 #Neglect the cross-section bending inertia
    XBINPUT.BeamMass[5,5] = 0.0 #Neglect the cross-section bending inertia
    XBINPUT.ForceStatic[-1,2] = 6e+05
    XBINPUT.ForceStatic[-1,1] = 6e+05
    XBINPUT.ForceStatic[-1,0] = 6e+03
    
    PosDefor, PsiDefor = NonlinearStatic.Solve_Py(XBINPUT, XBOPTS)
    
    "Initialise PyFSI variables"
    M = 10
    LeadingEdge = np.array([0.0, 0.5, 0.0])
    TrailingEdge = np.array([0.0, -0.5, 0.0])
    Section = FlatPlateReg(M,LeadingEdge,TrailingEdge)
    
    VelA_A = np.zeros(3)
    OmegaA_A = np.zeros(3)
    PosDotDef = np.zeros_like(PosDefor)
    PsiDotDef = np.zeros_like(PsiDefor)
    
    AeroGrid, AeroVels = Beam2UVLM.CoincidentGrid(PosDefor, PsiDefor, Section,\
                   VelA_A, OmegaA_A, PosDotDef, PsiDotDef,
                   XBINPUT)
    
    GridFile = 'AeroDef.dat' 
    
    fp = open(GridFile,'w')
    fp.write('VARIABLES= "X", "Y", "Z"\n')
    fp.write('ZONE I=%s, J=%s, DATAPACKING=BLOCK\n'\
              % (AeroGrid.shape[0],AeroGrid.shape[1]))
    for Var in range(3):
        for j in range(AeroGrid.shape[1]):
            for i in range(AeroGrid.shape[0]):
                fp.write('%f\t' % (AeroGrid[i,j,Var]))
                """if j == AeroGrid.shape[1]-1:
                    fp.write('\n')"""
            #END for j
        #END for i
        fp.write('\n')
    #END for Var
    
    
    
    