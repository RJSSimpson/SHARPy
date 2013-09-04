''' PyFSI.Utils.SectionInit
@brief      Functions to initialise wing sections for use with PyFSI.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       21/01/2013
@pre        None
@warning    None
'''

import numpy as np
import ctypes as ct

def FlatPlateReg(M,LeadingEdge,TrailingEdge):
    """@brief Regular interpolation between Triads LE and TE.
    @param M Number of chordwise panels.
    @param LeadingEdge Triad of LE point.
    @param TrailingEdge Triad of TE point.
    
    @details Sectional coordinates, i.e LE and TE, are be specified in the
    beam B-frame."""
    
    "check that x-component of LE and TE is zero"
    assert LeadingEdge[0] == 0.0, 'Sectional x-coord not zero!'
    assert TrailingEdge[0] == 0.0, 'Sectional x-coord not zero!'
    
    "Initialise empty array"
    Section = np.zeros((M+1,3), ct.c_double, 'C')
    
    "Get delta of panel displacements from LE to TE"
    Delta = (TrailingEdge - LeadingEdge) / float(M)
    
    for jSection in range(M+1):
        Section[jSection,:] = LeadingEdge + jSection * Delta
    
    return Section


if __name__ == '__main__':
    M = 5
    LeadingEdge = np.array([0.0, 0.5, 0.0])
    TrailingEdge = np.array([0.0, -0.5, 0.0])
    Section = FlatPlateReg(M, LeadingEdge, TrailingEdge)
    print(Section)
    pass