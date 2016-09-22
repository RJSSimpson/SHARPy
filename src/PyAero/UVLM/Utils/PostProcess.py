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

def GetCoeffs(VMOPTS, Forces, VMINPUT, VelA_G = None):
    """@brief Calculate force coefficients from UVLM solution.
    
    @returns force coefficients in frame defined by free-stream AoA."""
    
    # Declare temp traids.
    Psi = np.zeros((3))
    Coeff = np.zeros((3))
    
    # Sum all forces.
    for i in range(VMOPTS.M.value+1):
        for j in range(VMOPTS.N.value+1):
            Coeff[:] += Forces[i,j,:]
        #end for j
    #end for i
    
    # divide by denom.
    if VelA_G is not None:
        denom = 0.5*(np.linalg.norm(VMINPUT.U_infty-VelA_G))**2.0*VMINPUT.area
    else:
        denom = 0.5*np.linalg.norm(VMINPUT.U_infty)**2.0*VMINPUT.area
    
    Coeff[:] = Coeff[:]/denom
    
    # account for rotation of aerodynamic FoR (freestream velocity)
    Psi[0] = VMINPUT.alpha
    
    # get transformation matrix.
    CalphaG = Psi2TransMat(Psi)
    
    return np.dot(CalphaG,Coeff)



def WriteAeroTecHeader(FileName='AeroGrid.dat', Title='Default',
                          Variables=['X', 'Y', 'Z']):
    """@brief Opens Filename and writes header for tecplot ascii format.
    @return Open file object for writing to."""
    
    # Write title.
    FileObject = open(FileName,'w')
    FileObject.write('TITLE="%s"\n' % (Title))
    
    # If Variables are less than 3, write as grid file.
    if len(Variables) < 4:
        FileObject.write('FILETYPE=GRID ')
    elif len(Variables) >= 4:
        FileObject.write('FILETYPE=FULL ')
    
    # Write variables to the file.
    FileObject.write('VARIABLES=')
    for Var in range(len(Variables)):
        FileObject.write('"%s" ' % (Variables[Var]))
        if Var == len(Variables):
            FileObject.write('\n')
    #END for Var
    FileObject.write('\n')
    
    return FileObject


def WriteAeroTecZone(FileObject,
                       Name,
                       Zeta,
                       TimeStep=-1,
                       NumTimeSteps=0,
                       Time=0.0,
                       Variables = ['X', 'Y', 'Z'],
                       Text=True, 
                       Gamma = 0.0):
    """@brief Writes aero grid data to ascii file in tecplot ij format.
    @param FileObject Must be an open file object.
    @param Name Name of the block of data for teclplot [Surface/Wake]
    @param Zeta Lattice of points making up the block of data (M+1,N+1,3).
    @param Timestep -1 for static solution and >-1 thereafter.
    @param NumTimeSteps Total number of timesteps.
    @param Time Time at current step.
    @param Variables Pythonic list of variable names to plot, e.g
            ['X','Y','Z','Gamma'].
    @param Text Flag to write current time into output.
    @param Gamma Corresponding circulation strengths to plot (M,N).
    @details Tecplot ascii reader has a per-line character limit before which
    it looks for a new line."""
    
    STRANDID = 1
    if Name == 'Wake':
        STRANDID = 2
    
    FileObject.write('ZONE I=%d J=%d DATAPACKING=BLOCK T="%s: Timestep %d of %d"\nSOLUTIONTIME=%f STRANDID=%d\n' % \
                     (Zeta.shape[0],Zeta.shape[1],Name,TimeStep,\
                      NumTimeSteps,Time,STRANDID))
    if len(Variables) == 4: 
        FileObject.write('VARLOCATION = (4=CELLCENTERED)\n')
    
    for Var in range(len(Variables)):
        charcounter = 0
        for j in range(Zeta.shape[1]):
            for i in range(Zeta.shape[0]):
                if Var < 3:
                    FileObject.write('%f\t' % (Zeta[i,j,Var]))
                    charcounter += 9
                if (Var==3 and j<Zeta.shape[1]-1 and i<Zeta.shape[0]-1):
                    FileObject.write('%f\t' % (Gamma[i,j]))
                    charcounter += 9
                if charcounter > 900:
                    FileObject.write('\n')
                    charcounter = 0
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
    
def WriteUVLMtoTec(FileObject,
                     Zeta,
                     ZetaStar,
                     Gamma,
                     GammaStar,
                     TimeStep,
                     NumTimeSteps,
                     Time,
                     Text = False):
    """@brief Write surface and wake to open tecplot file object.
    @param FileObject Must be an open file object.
    @param Zeta Lattice of points on surface (M+1,N+1,3).
    @param ZetaStar Lattice of points wake (Mstar+1,N+1,3).
    @param Gamma Surface circulation strengths (M,N).
    @param GammaStar Wake circulation strengths (Mstar,N).
    @param Timestep -1 for static solution and >-1 thereafter.
    @param NumTimeSteps Total number of timesteps.
    @param Time Time at current step.
    @param Text Flag to write current time into output.
    @details If Text is true then the time only needs to be written once,
              and is done so when plotting the 'surface'.
    """
    Variables = ['X','Y','Z','Gamma']
    WriteAeroTecZone(FileObject,
                     'Surface',
                     Zeta,
                     TimeStep,
                     NumTimeSteps,
                     Time,
                     Variables,
                     Text,
                     Gamma)
    WriteAeroTecZone(FileObject,
                     'Wake',
                     ZetaStar,
                     TimeStep,
                     NumTimeSteps,
                     Time,
                     Variables,
                     False,
                     GammaStar)

if __name__ == '__main__':
    pass