'''@package PyAero.UVLM.Utils.PostProcess
@brief      Post-processing for UVLM.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       12/02/2013
@pre        None
@warning    None
'''

import numpy as np
from PyBeam.Utils.XbeamLib import Psi2TransMat

def GetCoeffs(VMOPTS, Forces, VMINPUT):
    """@brief Calculate force coefficients from UVLM solution."""
    
    "declare temp traids"
    Psi = np.zeros((3))
    Coeff = np.zeros((3))
    
    "sum all forces"
    for i in range(VMOPTS.M.value+1):
        for j in range(VMOPTS.N.value+1):
            Coeff[:] += Forces[i,j,:]
        #end for j
    #end for i
    
    "divide by denom"
    if VMINPUT.ZetaDotTest == 0.0:
        denom = 0.5*np.linalg.norm(VMINPUT.U_infty)**2.0*VMINPUT.area
    elif  VMINPUT.ZetaDotTest != 0.0:
        denom = 0.5*(np.linalg.norm(VMINPUT.U_infty)+VMINPUT.ZetaDotTest)**2.0*VMINPUT.area
    
    Coeff[:] = Coeff[:]/denom
    
    "account for rotation of aerodynamic FoR (freestream velocity)"
    Psi[0] = VMINPUT.alpha
    
    "get transformation matrix"
    CalphaG = Psi2TransMat(Psi)
    
    return np.dot(CalphaG,Coeff)



def WriteAeroTecHeader(FileName='AeroGrid.dat', Title='Default',\
                       Variables=['X', 'Y', 'Z']):
    """@brief Opens Filename and writes header for tecplot ascii format.
    @return open file object for writing to."""
    
    "write title"
    FileObject = open(FileName,'w')
    FileObject.write('TITLE="%s"\n' % (Title))
    
    "if Variables is less than 3, write as grid"
    if len(Variables) < 4:
        FileObject.write('FILETYPE=GRID ')
    elif len(Variables) >= 4:
        FileObject.write('FILETYPE=FULL ')
    
    "write variables"
    FileObject.write('VARIABLES=')
    for Var in range(len(Variables)):
        FileObject.write('"%s" ' % (Variables[Var]))
        if Var == len(Variables):
            FileObject.write('\n')
    #END for Var
    FileObject.write('\n')
    
    return FileObject


def WriteAeroTecZone(FileObject, Name, AeroGrid, TimeStep=-1, NumTimeSteps=0,\
                     Time=0.0, Variables=['X', 'Y', 'Z'], Text=True, Gamma = 0.0):
    """@brief Writes aero grid data to ascii file in tecplot ij format.
    @param FileObject Must be an open file object.
    @param Timestep -1 for static solution and >-1 thereafter."""
    
    STRANDID = 1 
    if Name == 'Wake':
        STRANDID = 2
    
    FileObject.write('ZONE I=%d J=%d DATAPACKING=BLOCK T="%s: Timestep %d of %d"\nSOLUTIONTIME=%f STRANDID=%d\n' % \
                     (AeroGrid.shape[0],AeroGrid.shape[1],Name,TimeStep+1,\
                      NumTimeSteps,Time,STRANDID))
    if len(Variables) == 4: 
        FileObject.write('VARLOCATION = (4=CELLCENTERED)\n')
    
    for Var in range(len(Variables)):
        for j in range(AeroGrid.shape[1]):
            for i in range(AeroGrid.shape[0]):
                if Var < 3:
                    FileObject.write('%f\t' % (AeroGrid[i,j,Var]))
                if (Var==3 and j<AeroGrid.shape[1]-1 and i<AeroGrid.shape[0]-1):
                    FileObject.write('%f\t' % (Gamma[i,j]))
            #END for j
        #END for i
        FileObject.write('\n')
    #END for Var
    
    if Text==True and TimeStep > -1:
        FileObject.write('TEXT T="time = %fs" CS=GRID3D AN=Headleft S=LOCAL ZN=%d X=0.0 Y=0.0 Z=-1.0\n' % (Time,TimeStep+1))
    
    return FileObject
    

def CloseAeroTecFile(FileObject):
    """TODO: could whole file writing process be encapsulated in a class?"""
    FileObject.close()

if __name__ == '__main__':
    pass