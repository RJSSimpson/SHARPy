'''@package PyAero.UVLM.UVLM
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
from DerivedTypesAero import VMopts, VMinput
import ctypes as ct
from XbeamLib import Psi2TransMat
import DerivedTypes
import PostProcess
import SharPySettings as Settings
import os
from VLM import InitSteadyExternalVels

Settings.OutputDir = os.getcwd() + '/'
Settings.OutputFileRoot = ''

def Run_Cpp_Solver_UVLM(VMOPTS,VMINPUT):
    """@brief UVLM solver with prescribed inputs."""
    
    "define reference line as in PyBeam"
    
    
    "init grid based section + definition of reference line"
    
    "init wake for unsteady solution"
    
    "init external velocities"  
    Uext = InitSteadyExternalVels(VMOPTS,VMINPUT)
    
    
    "Solve"
    UVLMLib.Cpp_Solver_VLM(Zeta, ZetaDot, Uext, ZetaStar, VMOPTS, \
                           Forces, Gamma, GammaStar)
    
    
    "Print tecplot file to check wake and grid etc"
    Variables=['X', 'Y', 'Z','Gamma']
    
    "write header"
    
    Filename = Settings.OutputDir + Settings.OutputFileRoot + 'AeroGrid.dat'
    
    FileObject = PostProcess.WriteAeroTecHeader(Filename,\
                                                'Default',\
                                                Variables)
    
    "write surface zone data"
    PostProcess.WriteAeroTecZone(FileObject, 'Surface', Zeta,\
                                -1, 0,\
                                0.0, Variables, False, Gamma)
    
    "write wake data"
    PostProcess.WriteAeroTecZone(FileObject, 'Wake', ZetaStar,\
                                -1, 0,\
                                0.0, Variables, False, GammaStar)
    
    "close file"
    PostProcess.CloseAeroTecFile(FileObject)
    
    "post process to get coefficients"
    return PostProcess.GetCoeffs(VMOPTS, Forces, VMINPUT)

if __name__ == '__main__':
    