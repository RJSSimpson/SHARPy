'''@package PyBeam.Utils.Input
@brief      Input functions to populate DerivedTypes for PyBeam simulations.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       28/11/2012
@pre        None
@warning    None
'''

import sys
import DerivedTypes

def Setup(XBINPUT, XBOPTS):
    """@brief Initial functionality adapted from FxInput_PyFx.f90."""
    
    "Check number of nodes per element and Gauss points"
    if XBINPUT.NumNodesElem != 2 and XBINPUT.NumNodesElem != 3:
        sys.stderr.write('Invalid number of nodes per element\n')
    elif XBINPUT.NumNodesElem == 2 and XBOPTS.NumGauss.value != 1:
        sys.stderr.write('Number of nodes per element (2) inconsistent with' +\
                         ' number of gauss points (%d)'\
                         %(XBOPTS.NumGauss.value) +\
                         ' Gauss points set to (1)\n')
        XBOPTS.NumGauss.value = 1
    elif XBINPUT.NumNodesElem == 3 and XBOPTS.NumGauss.value != 2:
        sys.stderr.write('Number of nodes per element (3) inconsistent with' +\
                         ' number of gauss points(%d)'\
                         %(XBOPTS.NumGauss.value)  +\
                         ' Gauss points set to (2)\n')
        XBOPTS.NumGauss.value = 2    
    sys.stderr.flush()
    
    "Apply the stiffness paramater to the stiffness matrix"
    XBINPUT.BeamStiffness = XBINPUT.BeamStiffness*XBINPUT.Sigma
        
    return XBINPUT, XBOPTS

def Elem(XBINPUT, XBOPTS, XBELEM):
    """@brief Set up Xbelem type. Initial functionality adapted from 
    FxInput_PyFx.f90."""
    
    "Check number of nodes per element and set default connectivity."
    assert (XBINPUT.NumNodesElem == 2 or XBINPUT.NumNodesElem == 3),\
            'Invalid number of nodes per element'
    if XBINPUT.NumNodesElem == 2:
        #TODO: check mapping of connectivity array from py to fortran
        pass
    elif XBINPUT.NumNodesElem == 3:
        #TODO: check mapping of connectivity array from py to fortran
        pass
    
    #TODO: calculate beam inverse stiffness
    #TODO: populate XBELEM. Stiff Mass InvStiff
    #TODO: Element orientation
    #TODO: member number (0 anyway)
    #DOING: test svn
    
    
if __name__ == '__main__':
    XBINPUT = DerivedTypes.Xbinput()
    XBINPUT.NumNodesElem = 3
    XBOPTS = DerivedTypes.Xbopts()
    XBINPUT, XBOPTS = Setup(XBINPUT, XBOPTS)
    XBINPUT, XBOPTS = Setup(XBINPUT, XBOPTS)