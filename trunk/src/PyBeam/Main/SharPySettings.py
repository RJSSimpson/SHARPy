'''@package PyBeam.Main.SharPySettings
@brief      System settings for SharPy
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       22/11/2012
@pre        None
@warning    None
'''

"""Directories"""

SharPyProjectDir = '/home/rjs10/SharPyProject/'

BeamLibDir = SharPyProjectDir + 'BeamLib/bin/'
BeamLibName = './BeamLib.so'

OutputDir =  SharPyProjectDir + 'Output/Temp/'
OutputFileRoot = 'Foo'

"""Structural Code Constants"""

MaxElNod = 3
DimMat = 18
RotEpsilon = 0.001 #rotations below this are linearised

if __name__ == '__main__':
    pass