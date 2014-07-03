"""@package PyBeam.Utils.BeamIO
@brief      Input/Output routines.
@author     Rob Simpson
@author     A. Da-Ronch (Liverpool)
@contact    r.simpson11@imperial.ac.uk
@version    1.0
@date       01/11/2012
@pre        None
@warning    None
"""

import sys, os
import ctypes as ct
import numpy as np
import SharPySettings as Settings
from PyBeam.Utils.Misc import iNode2iElem
import PyBeam.Utils.XbeamLib as XbeamLib
from PyBeam.Utils.BeamLib import Cbeam3_strainz

BeamPath = Settings.BeamLibDir + Settings.BeamLibName
BeamLib = ct.cdll.LoadLibrary(BeamPath)

f_fem_glob2loc_extract = BeamLib.__test_MOD_wrap_fem_glob2loc_extract

f_fem_glob2loc_extract.restype = None

def localStrains(PosDef, PsiDef,
                   PosIni, PsiIni,
                   XBELEM, NodeList,
                   SO3 = False):
    """@brief Approximate strains at the midpoint of nodes.
    
    @param PosDef Deformed nodal coordinates.
    @param PsiDef Deformed nodal rotation vectors.
    @param PosIni Initial (undeformed) nodal coordinates.
    @param PsiIni Initial (undeformed) nodal rotation vectors.
    @param XBELEM Xbeam element derived type containing element data.
    @param NodeList List of node numbers for strain approximation at adjacent
            midpoint in the direction of increasing node index,
            starts from zero:
            0---x---1---x---2---x---3-- ... --(NumNodes-1)
            Strains are calculated at the x to the right of selected nodes.
    @param SO3 Flag for use of the SO(3) manifold in CRV interpolation.
    
    @details Assumes that all nodes are either two- or three- noded. 
    @warning Untested.
    """
    
    strains = np.zeros((len(NodeList),6))
    NumNodes = PosDef.shape[0]-1
    NumNodesElem = XBELEM.NumNodes[0]
    iOut = 0 # Output index
    for iNode in NodeList:
        if iNode == NumNodes or iNode == -1:
            raise ValueError("Midpoint strain requested beyond final node.")
        
        iElem, iiElem = iNode2iElem(iNode, NumNodes, NumNodesElem)
        
        if NumNodesElem == 2:
            if SO3 == True:
                # Consistent calculation of midpoint strains using SO(3) manifold.
                raise NotImplementedError("SO(3) manifold strains.")
            else:
                # Using fortran routines
                R0 = np.zeros((2,6))
                R0[0,:3] = PosIni[iNode,:]
                R0[0,3:] = PsiIni[iElem,iiElem,:]
                R0[1,:3] = PosIni[iNode+1,:]
                R0[1,3:] = PsiIni[iElem,iiElem+1,:]
                
                Ri = np.zeros((2,6))
                Ri[0,:3] = PosDef[iNode,:]
                Ri[0,3:] = PsiDef[iElem,iiElem,:]
                Ri[1,:3] = PosDef[iNode+1,:]
                Ri[1,3:] = PsiDef[iElem,iiElem+1,:]
            
                # Midpoint
                z = 0.0
                
                strainz = Cbeam3_strainz(R0, Ri, z)
            # END if SO3            
            
        elif NumNodesElem == 3:
            if SO3 == True:
                raise NotImplementedError("Only available for 2-noded just now.")
            else:
                # Using fortran routines
                R0 = np.zeros((3,6))
                R0[0,:3] = PosIni[2*iElem,:]
                R0[1,:3] = PosIni[2*iElem+2,:]
                R0[2,:3] = PosIni[2*iElem+1,:]
                R0[:,3:] = PsiIni[iElem]
                
                Ri = np.zeros((3,6))
                Ri[0,:3] = PosDef[2*iElem,:]
                Ri[1,:3] = PosDef[2*iElem+2,:]
                Ri[2,:3] = PosDef[2*iElem+1,:]
                Ri[:,3:] = PsiDef[iElem]
                
                if iiElem == 0:
                    z = -0.5 # mid-point adjacent to node iNode
                elif iiElem == 1:
                    z = 0.5 # mid-point adjacent to node iNode
                elif iiElem == 2:
                    raise ValueError("iiElem only 0, 1 for adjacent midpoint" +
                                     "calculation.")
                else:
                    raise ValueError("iiElem can only take 0, 1 or 2 for " +
                                     "3-noded elements.")
                
                strainz = Cbeam3_strainz(R0, Ri, z)
            # END if SO(3)
        else:
            raise ValueError("Only 2- or 3-noded elements are supported.")
        
        strains[iOut,:] = strainz
        
        iOut += 1
    # END for iElem
    return strains
    
def localElasticForces(PosDef, PsiDef, 
                          PosIni, PsiIni,
                          XBELEM, NodeList):
    """@brief Approximate shear force and moments from local strain and
    stiffness matrix.
    
    @param PosDef Deformed nodal coordinates.
    @param PsiDef Deformed nodal rotation vectors.
    @param XBELEM Xbeam element derived type containing element data.
    @param NodeList List of element numbers for strain approximation,
            starts from zero.
    @warning Untested.
    """
    
    F = np.zeros((len(NodeList),6))
    strains = localStrains(PosDef, PsiDef, 
                           PosIni, PsiIni,
                           XBELEM, NodeList)
    iOut = 0 # Output index
    for iNode in NodeList:
        iElem, iiElem = iNode2iElem(iNode, PosDef.shape[0]-1, XBELEM.NumNodes[0])
        del iiElem
        elemStrain = strains[iOut,:]
        elemK = XBELEM.Stiff[iElem*6:(iElem+1)*6,:]
        F[iOut,:] = np.dot(elemK,elemStrain)
        iOut += 1
    # END of iElem
    return F

#-------------------------------------------------------------------------------
# START - DisplayProgressCounter
#-------------------------------------------------------------------------------
def DisplayProgressCounter(ni,nt):
    """
    """
    sys.stdout.write('\r    %3d/%3d' %(ni,nt)); sys.stdout.flush()
    
    return ni+1
#-------------------------------------------------------------------------------
# END - DisplayProgressCounter
#-------------------------------------------------------------------------------

def OutputElems(NumElems,TotNumNodes,XBELEM,PosDefor,PsiDefor,ofile,WriteMode):
    """@brief Writes beam displacements and rotations to an output file.
    
    @param NumElems Number of elements in the model.
    @param TotNumNodes Number of nodes, pythonic int.
    @param XBELEM Xbelem derived type.
    @param PosDefor Array of nodal position vectors.
    @param PsiDefor Array of nodal orientation CRVs.
    @param ofile Pythonic file object to write to.
    @param WriteMode Write mode for open() function"""
    
    PosGlob     = np.zeros((Settings.MaxElNod*NumElems,3),\
                              dtype=np.double,\
                              order="Fortran")
    NumNE_array = (ct.c_int *NumElems)()
    
    f_fem_glob2loc_extract( \
        ct.byref(ct.c_int(NumElems)), \
        ct.byref(ct.c_int(TotNumNodes)), \
        XBELEM.NumNodes.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBELEM.MemNo.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBELEM.Conn.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBELEM.Master.ctypes.data_as(ct.POINTER(ct.c_int)), \
        XBELEM.Length.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.PreCurv.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.Psi.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.Vector.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.Mass.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.Stiff.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.InvStiff.ctypes.data_as(ct.POINTER(ct.c_double)), \
        XBELEM.RBMass.ctypes.data_as(ct.POINTER(ct.c_double)),\
        PosDefor.ctypes.data_as(ct.POINTER(ct.c_double)), \
        PsiDefor.ctypes.data_as(ct.POINTER(ct.c_double)), \
        PosGlob.ctypes.data_as(ct.POINTER(ct.c_double)), \
        ct.byref(NumNE_array) )
    
    fp = open(ofile,WriteMode)
    i1 = 0
    for iElem in range(NumElems):
        for i in range(NumNE_array[iElem]):
            fp.write('%d %d %13.6e %13.6e %13.6e %13.6e %13.6e %13.6e\n'
                %(iElem+1,i+1,PosGlob[i1,0],PosGlob[i1,1],PosGlob[i1,2],
                  PsiDefor[iElem,i,0],PsiDefor[iElem,i,1],PsiDefor[iElem,i,2]))
            i1 += 1
    fp.close()
    
    return

def WriteUndefGeometry(NumElems,TotNumNodes,XBELEM,
                          PosIni,PsiIni,RootFile,WriteMode = 'a'):
    """@brief Create *_UndefGeom.dat file and write output of 
    fem_glob2loc_extract from lib_fem.f90.
    
    Default behaviour is to append to existing file
    TotNumNodes is a Python int equal to NumNodes_tot.value"""

    ofile = RootFile + '_und.dat'
    fp = open(ofile,WriteMode)
    fp.write('TITLE="Undeformed initial geometry"\n')
    fp.write('VARIABLES="iElem" "iNode" "Px" "Py" "Pz" "Rx" "Ry" "Rz"\n')
    fp.close()

    OutputElems(NumElems,TotNumNodes,XBELEM,
                PosIni,PsiIni,ofile,WriteMode)

    return

#-------------------------------------------------------------------------------
# START - Write_force_File
#-------------------------------------------------------------------------------
def Write_force_File(fp,Time,ForceTime,ForcedVel,ForcedVelDot):
    """
    """

    nrl = np.size(ForcedVel,0)
    
    fp.write('TITLE="Forcing velocity/acceleration of support"\n')
    fp.write('VARIABLES="T" "Ft" "V1" "V2" "V3" "V4" "V5" "V6" "A1" "A2" "A3" "A4" "A5" "A6"\n')
    fp.write('ZONE I=%d T=""\n' %(nrl))
    
    for i1 in range(nrl):
        fp.write(' %12.5e %12.5e ' %(Time[i1],ForceTime[i1]))
        for i2 in range(6):
            fp.write('%12.5e ' %(ForcedVel[i1,i2]))
        for i2 in range(6):
            fp.write('%12.5e ' %(ForcedVelDot[i1,i2]))
        fp.write('\n')
    
    return
#-------------------------------------------------------------------------------
# END - Write_force_File
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# START - Write_dyn_File
#-------------------------------------------------------------------------------
def Write_dyn_File(fp,Time,PosPsiTime):
    """ @brief Writes _dyn file based on data from PosPsiTime.
    """
    
    nr1 = np.size(PosPsiTime,0); nr2 = np.size(PosPsiTime,1)

#     fp.write('TITLE="Wing-tip dynamic response"\n')
#     fp.write('VARIABLES="T" "PX" "PY" "PZ" "RX" "RY" "RZ"\n')
#     fp.write('ZONE I=%d T="%s"\n' %(nr1,'dyn'))

    for i1 in range(nr1):
        fp.write(' %12.5e ' %(Time[i1]))
        for i2 in range(nr2):
            fp.write('%12.5e ' %(PosPsiTime[i1,i2]))
        fp.write('\n')
    
    return
#-------------------------------------------------------------------------------
# END - Write_dyn_File
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# START - Write_vel_File
#-------------------------------------------------------------------------------
def Write_vel_File(fp,VelocTime):
    """
    """

    nr1 = np.size(VelocTime,0); nr2 = np.size(VelocTime,1)

    for i1 in range(nr1):
        for i2 in range(nr2):
            fp.write(' %12d %12d %12.5e\n' %(i1+1, i2+1, VelocTime[i1,i2]))

    return
#-------------------------------------------------------------------------------
# END - Write_vel_File
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# START - Write_shape_File
#-------------------------------------------------------------------------------
def Write_shape_File(fp,TotNrSteps,TotNrNodes,Time,DynOut):
    """
    """
    
#     fp.write('TITLE="Deformed shape function of time"\n')
#     fp.write('VARIABLES="T" "PX" "PY" "PZ"\n')
    
    for i1 in range(TotNrSteps):
#         fp.write('ZONE I=%d T="iStep: %d"\n' %(TotNrNodes,i1))
        for i2 in range(TotNrNodes):
            i3 = i1*TotNrNodes + i2
            fp.write(' %12.5e %12.5e %12.5e %12.5e\n' \
                %(Time[i1], DynOut[i3,0], DynOut[i3,1], DynOut[i3,2]))
    
    return
#-------------------------------------------------------------------------------
# END - Write_shape_File
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# START - Write_rigid_File
#-------------------------------------------------------------------------------
def Write_rigid_File(fp,Time,RefVel,RefVelDot):
    """
    """
    
    nr1 = np.size(RefVel,0); nr2 = np.size(RefVel,1)
    
#     fp.write('TITLE="Rigid body velocity/acceleration"\n')
#     fp.write('VARIABLES="T" "V1" "V2" "V3" "V4" "V5" "V6" "A1" "A2" "A3" "A4" "A5" "A6"\n')
#     fp.write('ZONE I=%d T=""\n' %(nr1))

    for i1 in range(nr1):
        fp.write(' %12.5e ' %(Time[i1]))
        for i2 in range(nr2):
            fp.write('%12.5e ' %(RefVel[i1,i2]))
        for i2 in range(nr2):
            fp.write('%12.5e ' %(RefVelDot[i1,i2]))
        fp.write('\n')
    
    return
#-------------------------------------------------------------------------------
# END - Write_rigid_File
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# START - Write_ComplexVector
#-------------------------------------------------------------------------------
def Write_ComplexVector(fp,Vec):
    """
    """
    
    nr1 = np.size(Vec)
    
    for i1 in range(nr1):
        fp.write(' %12.5e %12.5e\n' %(Vec[i1].real,Vec[i1].imag))
    
    return
#-------------------------------------------------------------------------------
# END - Write_ComplexVector
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# START - Write_Eigsol
#-------------------------------------------------------------------------------
def Write_Eigsol(fp,Vec,Mat):
    """
    """
    
    nr1 = np.size(Mat,1); nr2 = np.size(Mat,0)
    
    fp.write('# Eigenfrequencies\n')
    for i1 in range(nr1):
        fp.write(' %12.5e %12.5e |' %(np.real(Vec[i1]),np.imag(Vec[i1])))
    fp.write('\n# Eigenmodes\n')
    for i2 in range(nr2):
        for i1 in range(nr1):
            fp.write(' %12.5e %12.5e |' %(np.real(Mat[i2,i1]),\
                                          np.imag(Mat[i2,i1])))
        fp.write('\n')
    
    return
#-------------------------------------------------------------------------------
# END - Write_Eigsol
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# START - dump_partitioned_residual
#-------------------------------------------------------------------------------
def dump_partitioned_residual(res,n1,n2):
    """
    res: residual partitioned in structure+fluid, dim(res)=n1+n2
    n1 : number of structure dofs
    n2 : number of fluid dofs
    """

    tol = 1.e-10; tol1 = 1.e-03

    for i1 in range(n1+n2):
        if(i1 == 0):
            sys.stdout.write('DEBUG:    Structural system residual (zeroed if .lt. %e)\n' %(tol))
        elif(i1 == n1):
            sys.stdout.write('DEBUG:    Fluid system residual (zeroed if .lt. %e)\n' %(tol))
        if(abs(res[i1]) < tol):
            sys.stdout.write('DEBUG:     %12.5e <--- (.lt. %e)\n' %(0.,tol))
        else:
            sys.stdout.write('DEBUG:     %12.5e' %(res[i1]))
            if(abs(res[i1]) > tol1):
                sys.stdout.write(' <+++ (.gt. %e)' %(tol1))
            sys.stdout.write('\n')

    return
#-------------------------------------------------------------------------------
# END - dump_partitioned_residual
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# START - Read_nastran_file
#-------------------------------------------------------------------------------
def Read_nastran_file(NastranFile):
    """
    """
    
    sys.stdout.write('\n *** Reading input file: %s\n' %(NastranFile))

    FIELD = 8
    
    # Output structure
    OutData = {}

    NLPARAM = {'ID': [],
               'NINC': [],
               'DT': [],
               'KMETHOD': [],
               'KSTEP': [],
               'MAXITER': [],
               'CONV': [],
               'INTOUT': [],
               'EPSU': [],
               'EPSP': [],
               'EPSW': [],
               'MAXDIV': [],
               'MAXQN': [],
               'MAXLS': [],
               'FSTRESS': [],
               'LSTOL': [],
               'MAXBIS': [],
               'MAXR': [],
               'RTOLB': []
              }

    AERO = {'ACSID': [],
            'VELOCITY': [],
            'REFC': [],
            'RHOREF': [],
            'SYMXZ': [],
            'SYMXY': []
           }

    GUST = {'SID': [],
            'DLOAD': [],
            'WG': [],
            'X0': [],
            'V': []
           }

    GRID = {'ID': [],
            'CP': [],
            'X1': [],
            'X2': [],
            'X3': [],
            'CD': [],
            'PS': [],
            'SEID': []
           }
    

    CBAR = {'EID': [],
            'PID': [],
            'GA': [],
            'GB': [],
            'X1': [],
            'X2': [],
            'X3': [],
            'OFFT': []
           }

    EIGR = {'SID': [],
            'METHOD': [],
            'F1': [],
            'F2': [],
            'NE': [],
            'ND': [],
            'NORM': [],
            'G': [],
            'C': []
           }

    TSTEP = {'SID': [],
             'N1': [],
             'DT1': [],
             'NO1': []
            }

    FORCE = {'FORCE': [],
             'SID': [],
             'G': [],
             'CID': [],
             'F': [],
             'N1': [],
             'N2': [],
             'N3': []
            }
    
    if(os.path.isfile(NastranFile) == False):
        sys.stderr.write('\n\n\tERROR - file %s not found!\n' %(NastranFile))
    
    fp = open(NastranFile,'r')
    line = fp.readline(); sys.stdout.write('%s' %(line))
    
    while(len(line) > 0):
        
        data = line.split(); INFO = 0
        
        # remove comment and empty lines
        if( (line[0] != '$') & (data != []) ):
            
            if(len(data[0])>FIELD):
                data = line.split(','); INFO = 1
            if(INFO==0): data = line
            
            CARD = field_parser(extract_field_nas(data,0,INFO),'s');
            
            #-------------------------------------------------------------------
            # NLPARAM card
            #-------------------------------------------------------------------
            if(CARD == 'NLPARAM'):
                
                NLPARAM['ID']      = field_parser(extract_field_nas(data,1,INFO),'d')
                NLPARAM['NINC']    = field_parser(extract_field_nas(data,2,INFO),'d')
                NLPARAM['DT']      = field_parser(extract_field_nas(data,3,INFO),'f')
                NLPARAM['KMETHOD'] = field_parser(extract_field_nas(data,4,INFO),'s')
                NLPARAM['KSTEP']   = field_parser(extract_field_nas(data,5,INFO),'d')
                NLPARAM['MAXITER'] = field_parser(extract_field_nas(data,6,INFO),'d')
                NLPARAM['CONV']    = field_parser(extract_field_nas(data,7,INFO),'s')
                NLPARAM['INTOUT']  = field_parser(extract_field_nas(data,8,INFO),'s')
                
                # try to read continuation
                last = field_parser(extract_field_nas(data,9,INFO),'s')
                
                if( len(last) > 0 ):
                    
                    line = fp.readline(); sys.stdout.write('%s' %(line))
                    
                    data = line.split(); INFO = 0
                    if(len(data[0])>FIELD):
                        data = line.split(','); INFO = 1
                    if(INFO==0): data = line
                    
                    first = field_parser(extract_field_nas(data,0,INFO),'s')
                    
                    if( first == last ):
                        
                        NLPARAM['EPSU']    = field_parser(extract_field_nas(data,1,INFO),'f')
                        NLPARAM['EPSP']    = field_parser(extract_field_nas(data,2,INFO),'f')
                        NLPARAM['EPSW']    = field_parser(extract_field_nas(data,3,INFO),'f')
                        NLPARAM['MAXDIV']  = field_parser(extract_field_nas(data,4,INFO),'d')
                        NLPARAM['MAXQN']   = field_parser(extract_field_nas(data,5,INFO),'d')
                        NLPARAM['MAXLS']   = field_parser(extract_field_nas(data,6,INFO),'d')
                        NLPARAM['FSTRESS'] = field_parser(extract_field_nas(data,7,INFO),'f')
                        NLPARAM['LSTOL']   = field_parser(extract_field_nas(data,8,INFO),'f')
                        
                # try to read continuation
                last = field_parser(extract_field_nas(data,9,INFO),'s')
                
                if( len(last) > 0 ):
                    
                    line = fp.readline(); sys.stdout.write('%s' %(line))
                    
                    data = line.split(); INFO = 0
                    if(len(data[0])>FIELD):
                        data = line.split(','); INFO = 1
                    if(INFO==0): data = line
                    
                    first = field_parser(extract_field_nas(data,0,INFO),'s')
                    
                    if( first == last ):
                        
                        NLPARAM['MAXBIS'] = field_parser(extract_field_nas(data,1,INFO),'d')
                        NLPARAM['MAXR']   = field_parser(extract_field_nas(data,5,INFO),'f')
                        NLPARAM['RTOLB']  = field_parser(extract_field_nas(data,7,INFO),'f')
                
                OutData['NLPARAM'] = NLPARAM
                line = fp.readline(); sys.stdout.write('%s' %(line))

            #-------------------------------------------------------------------
            # AERO card
            #-------------------------------------------------------------------
            elif(CARD == 'AERO'):
                
                AERO['ACSID']    = field_parser(extract_field_nas(data,1,INFO),'d')
                AERO['VELOCITY'] = field_parser(extract_field_nas(data,2,INFO),'f')
                AERO['REFC']     = field_parser(extract_field_nas(data,3,INFO),'f')
                AERO['RHOREF']   = field_parser(extract_field_nas(data,4,INFO),'f')
                AERO['SYMXZ']    = field_parser(extract_field_nas(data,5,INFO),'d')
                AERO['SYMXY']    = field_parser(extract_field_nas(data,6,INFO),'d')
                
                OutData['AERO'] = AERO
                line = fp.readline(); sys.stdout.write('%s' %(line))

            #-------------------------------------------------------------------
            # GUST card
            #-------------------------------------------------------------------
            elif(CARD == 'GUST'):
                
                GUST['SID']   = field_parser(extract_field_nas(data,1,INFO),'d')
                GUST['DLOAD'] = field_parser(extract_field_nas(data,2,INFO),'d')
                GUST['WG']    = field_parser(extract_field_nas(data,3,INFO),'f')
                GUST['X0']    = field_parser(extract_field_nas(data,4,INFO),'f')
                GUST['V']     = field_parser(extract_field_nas(data,5,INFO),'f')
                
                OutData['GUST'] = GUST
                line = fp.readline(); sys.stdout.write('%s' %(line))

            #-------------------------------------------------------------------
            # GRID card
            #-------------------------------------------------------------------
            elif(CARD == 'GRID'):
                
                GRID['ID'].append(field_parser(extract_field_nas(data,1,INFO),'d'))
                GRID['CP'].append(field_parser(extract_field_nas(data,2,INFO),'d'))
                GRID['X1'].append(field_parser(extract_field_nas(data,3,INFO),'f'))
                GRID['X2'].append(field_parser(extract_field_nas(data,4,INFO),'f'))
                GRID['X3'].append(field_parser(extract_field_nas(data,5,INFO),'f'))
                GRID['CD'].append(field_parser(extract_field_nas(data,6,INFO),'d'))
                GRID['PS'].append(field_parser(extract_field_nas(data,7,INFO),'d'))
                GRID['SEID'].append(field_parser(extract_field_nas(data,8,INFO),'d'))
                
                OutData['GRID'] = GRID
                line = fp.readline(); sys.stdout.write('%s' %(line))

            #--------------------------------------------------------------------------
            # CBAR card
            #--------------------------------------------------------------------------
            elif(CARD == 'CBAR'):
                
                CBAR['EID'].append( field_parser(extract_field_nas(data,1,INFO),'d'))
                CBAR['PID'].append( field_parser(extract_field_nas(data,2,INFO),'d'))
                CBAR['GA'].append(  field_parser(extract_field_nas(data,3,INFO),'d'))
                CBAR['GB'].append(  field_parser(extract_field_nas(data,4,INFO),'d'))
                CBAR['X1'].append(  field_parser(extract_field_nas(data,5,INFO),'f'))
                CBAR['X2'].append(  field_parser(extract_field_nas(data,6,INFO),'f'))
                CBAR['X3'].append(  field_parser(extract_field_nas(data,7,INFO),'f'))
                CBAR['OFFT'].append(field_parser(extract_field_nas(data,8,INFO),'s'))
                
                OutData['CBAR'] = CBAR
                line = fp.readline(); sys.stdout.write('%s' %(line))

            #--------------------------------------------------------------------------
            # EIGR card
            #--------------------------------------------------------------------------
            elif(CARD == 'EIGR'):
                
                EIGR['SID'].append(   field_parser(extract_field_nas(data,1,INFO),'d'))
                EIGR['METHOD'].append(field_parser(extract_field_nas(data,2,INFO),'s'))
                EIGR['F1'].append(    field_parser(extract_field_nas(data,3,INFO),'f'))
                EIGR['F2'].append(    field_parser(extract_field_nas(data,4,INFO),'f'))
                EIGR['NE'].append(    field_parser(extract_field_nas(data,5,INFO),'d'))
                EIGR['ND'].append(    field_parser(extract_field_nas(data,6,INFO),'d'))

                # try to read continuation
                last = field_parser(extract_field_nas(data,9,INFO),'s')
                
                if( len(last) > 0 ):
                    
                    line = fp.readline(); sys.stdout.write('%s' %(line))
                    
                    data = line.split(); INFO = 0
                    if(len(data[0])>FIELD):
                        data = line.split(','); INFO = 1
                    if(INFO==0): data = line
                    
                    first = field_parser(extract_field_nas(data,0,INFO),'s')

                    if( first == last ):
                        
                        EIGR['NORM'] = field_parser(extract_field_nas(data,1,INFO),'s')
                        EIGR['G']    = field_parser(extract_field_nas(data,2,INFO),'d')
                        EIGR['C']    = field_parser(extract_field_nas(data,2,INFO),'d')
                
                OutData['EIGR'] = EIGR
                line = fp.readline(); sys.stdout.write('%s' %(line))

            #--------------------------------------------------------------------------
            # TSTEP card
            #--------------------------------------------------------------------------
            elif(CARD == 'TSTEP'):
                
                TSTEP['SID'] =      field_parser(extract_field_nas(data,1,INFO),'d')
                TSTEP['N1'].append( field_parser(extract_field_nas(data,2,INFO),'d'))
                TSTEP['DT1'].append(field_parser(extract_field_nas(data,3,INFO),'f'))
                TSTEP['NO1'].append(field_parser(extract_field_nas(data,4,INFO),'d'))
                
                # try to read continuation
                iskip = 1
                while(iskip == 1):
                    
                    line = fp.readline(); sys.stdout.write('%s' %(line))
                    
                    # assume data are written in the same format as primary line
                    data = line.split()
                    if(INFO == 1): data = line.split(',')
                    else: data = line
                    
                    # check if first two fields are empty
                    if(not field_parser(extract_field_nas(data,0,INFO),'s') and \
                       not field_parser(extract_field_nas(data,1,INFO),'s') and \
                           field_parser(extract_field_nas(data,2,INFO),'d')):
                        
                        sys.stdout.write('%s' %(line))

                        TSTEP['N1'].append( field_parser(extract_field_nas(data,2,INFO),'d'))
                        TSTEP['DT1'].append(field_parser(extract_field_nas(data,3,INFO),'f'))
                        TSTEP['NO1'].append(field_parser(extract_field_nas(data,4,INFO),'d'))
                        
                    else:
                        iskip = 0
                
                OutData['TSTEP'] = TSTEP

            #--------------------------------------------------------------------------
            # FORCE card
            #--------------------------------------------------------------------------
            elif(CARD == 'FORCE'):
                
                FORCE['SID'].append(field_parser(extract_field_nas(data,1,INFO),'d'))
                FORCE['G'].append(  field_parser(extract_field_nas(data,2,INFO),'d'))
                FORCE['CID'].append(field_parser(extract_field_nas(data,3,INFO),'d'))
                FORCE['F'].append(  field_parser(extract_field_nas(data,4,INFO),'f'))
                FORCE['N1'].append( field_parser(extract_field_nas(data,5,INFO),'f'))
                FORCE['N2'].append( field_parser(extract_field_nas(data,6,INFO),'f'))
                FORCE['N3'].append( field_parser(extract_field_nas(data,7,INFO),'f'))
               
                OutData['FORCE'] = FORCE
                line = fp.readline(); sys.stdout.write('%s' %(line))

            else:
                
                # move to next line
                line = fp.readline(); sys.stdout.write('%s' %(line))
                
        else:
            
            # move to next line
            line = fp.readline(); sys.stdout.write('%s' %(line))
        
    fp.close()
    
    #---------------------------------------------------------------------------
    # Check consistency
    #---------------------------------------------------------------------------
    if( (GUST != {}) & (AERO != {}) ):
        if( (GUST['V']!=[]) & (AERO['VELOCITY']!=[]) ):
            if( GUST['V']!=AERO['VELOCITY'] ):
                sys.stderr.write('\nERROR: GUST.V=%f must be equal to AERO.VELOCITY=%f\n' 
                    %(GUST['V'],AERO['VELOCITY']))
                exit(-1)
    
    sys.stdout.write('\n *** Finished reading input file: %s\n' %(NastranFile))
    
    return OutData
#-------------------------------------------------------------------------------
# END - Read_nastran_file
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# START - field_parser
#-------------------------------------------------------------------------------
def field_parser(line,flag):
    """
    flag: 'd' --> integer
          'f' --> float
          's' --> string
    
    Parse and convert to specified type
    """
    
    num = []
    
    if(line == []): tmp = []
    else: tmp = line.split()
    
    if(len(tmp) > 0):
        if  (flag == 'd'): num = int(str(tmp[0]))
        elif(flag == 'f'):
            # can be largely improved and reduced in lines
            if( ('e' in tmp[0]) | ('E' in tmp[0]) ):
                num = float(str(tmp[0]))
            else:
                isign = -1
                for i in range(len(tmp[0])):
                    if( (i>0) & ((tmp[0][i]=='+') | (tmp[0][i]=='-')) ):
                        isign = i
                if(isign != -1):
                    num = float(str(tmp[0][:isign]))
                    num *= 10**( int(str(tmp[0][isign:])) )
                else:
                    num = float(str(tmp[0]))
        elif(flag == 's'): num = str(tmp[0])
        
    return num
#-------------------------------------------------------------------------------
# END - field_parser
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# START - extract_field_nas
#-------------------------------------------------------------------------------
def extract_field_nas(line,nr,flag):
    """
    """
    out = []
    
    if(flag==0):
        FIELD = 8
        out = line[FIELD*nr:FIELD*(nr+1)]
    else:
        if(nr<len(line)): out = line[nr]
    
    return out
#-------------------------------------------------------------------------------
# END - extract_field_nas
#-------------------------------------------------------------------------------

