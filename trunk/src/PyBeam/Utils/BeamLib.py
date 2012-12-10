'''@package PyBeam.Utils.BeamLib
@brief      Loads the f90 subroutines.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk 
@version    0.0
@date       10/12/2012
@pre        None
@warning    None
'''

import ctypes as ct #http://docs.python.org/3.2/library/ctypes.html
import SharPySettings as Settings

BeamPath = Settings.BeamLibDir + Settings.BeamLibName
BeamLib = ct.cdll.LoadLibrary(BeamPath)

f_input_setup               = BeamLib.__test_MOD_wrap_input_setup
f_input_elem                = BeamLib.__test_MOD_wrap_input_elem
f_input_node                = BeamLib.__test_MOD_wrap_input_node
f_xbeam_undef_geom          = BeamLib.__test_MOD_wrap_xbeam_undef_geom
f_xbeam_undef_dofs          = BeamLib.__test_MOD_wrap_xbeam_undef_dofs
f_cbeam3_solv_nlnstatic      = BeamLib.__test_MOD_wrap_cbeam3_solv_nlnstatic

"""ctypes does not check whether the correct number OR type of input arguments
are passed to each of these functions - great care must be taken to ensure the
number of arguments and the argument types are correct.
TODO: developer provided argtypes which must be a sequence of C data types.
http://docs.python.org/2/library/ctypes.html"""

f_input_setup.restype               = None
f_input_elem.restype                = None
f_input_node.restype                = None
f_xbeam_undef_geom.restype          = None
f_xbeam_undef_dofs.restype          = None
f_cbeam3_solv_nlnstatic.restype     = None

if __name__ == '__main__':
    pass