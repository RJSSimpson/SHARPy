'''@package PyBeam.Main.SharPySettings
@brief      System settings for SharPy
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       22/11/2012
@pre        None
@warning    None
'''
import getpass
import sys


# Directories.
userid = getpass.getuser()

SharPyProjectDir = '/home/' + userid + '/git/SHARPy/'

BeamLibDir = SharPyProjectDir + 'BeamLib/bin/'
BeamLibName = './BeamLib.so'

UVLMLibDir = SharPyProjectDir + 'UVLMLib/Debug/'
UVLMLibName = './UVLMLib.so'

OutputDir =  SharPyProjectDir + 'output/temp/'
OutputFileRoot = 'Foo'


# Python path.
sys.path.append(SharPyProjectDir + 'src')
sys.path.append(SharPyProjectDir + 'src/PyBeam')
sys.path.append(SharPyProjectDir + 'src/PyBeam/Utils')
sys.path.append(SharPyProjectDir + 'src/PyBeam/Solver')
sys.path.append(SharPyProjectDir + 'src/PyBeam/Main')
sys.path.append(SharPyProjectDir + 'src/PyFSI')
sys.path.append(SharPyProjectDir + 'src/PyFSI/Utils')
sys.path.append(SharPyProjectDir + 'src/PyAero')
sys.path.append(SharPyProjectDir + 'src/PyAero/UVLM')
sys.path.append(SharPyProjectDir + 'src/PyAero/UVLM/Utils')
sys.path.append(SharPyProjectDir + 'src/PyAero/UVLM/Solver')
sys.path.append(SharPyProjectDir + 'src/PyCoupled')
sys.path.append(SharPyProjectDir + 'src/PyCoupled/Utils')

# Structural Code Constants.
MaxElNod = 3
DimMat = 24 # Memory allocated for sparse matrix storage in fortran solver.
RotEpsilon = 0.001 # Rotations below this are linearised.


# Tecplot.
PlotTec = True

# Output options
WriteOut = True

if __name__ == '__main__':
    pass