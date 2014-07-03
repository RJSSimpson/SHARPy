'''
Created on 5 Nov 2013

@author: rjs10
'''

import ctypes as ct
from collections import OrderedDict
from PyBeam.Utils.BeamIO import localElasticForces, localStrains
import re

def OpenOutFile(writeDict, XBOPTS, Settings):
    """@brief Open an _out.dat file and write the header.
    
    @param writeDict An ordered dict of outputs to write.
    @param XBOPTS PyBeam options from which the solution code is obtained.
    @param Settings Current SharpySettings module. 
    @returns fp An open _out.dat file object. 
    """
    
    ofile = Settings.OutputDir + \
            Settings.OutputFileRoot + \
            '_SOL' + str(XBOPTS.Solution.value) + '_out.dat'
    fp = open(ofile,'w')
    fp.write("{:<14}".format("Time"))
    for output in writeDict.keys():
        fp.write("{:<14}".format(output))
    fp.write("\n")
    fp.flush()
    
    return fp


def WriteToOutFile(writeDict,
                     fp,
                     Time_iStep,
                     PosDefor,
                     PsiDefor,
                     PosIni,
                     PsiIni,
                     XBELEM,
                     ctrlSurf = None,
                     mpcCont = None):
    """@brief Write simulation outputs to file object.
    
    @param writeDict An ordered dict of outputs to write.
    @param fp An open _out.dat file object.
    @param Time_iStep Current time.
    @param PosDefor Current nodal displacements.
    @param PsiDefor Current element CRVs.
    @param PosIni Initial nodal displacements.
    @param PsiIni Initial element CRVs.
    @param XBELEM Element info derived type from PyBeam.
    @param ctrlSurf Control surface derived type.
    @param mpcCont MPC controller.
    """
    
    fp.write("{:<14,e}".format(Time_iStep))
    for myStr in writeDict.keys():
        if re.search(r'^R_[xyz]',myStr):
            if re.search(r'^R_[xyz]_[0-9]', myStr):
                index = int(myStr[4])
            elif re.search(r'root', myStr):
                index = 0
            elif re.search(r'tip', myStr):
                index = -1
            else:
                raise IOError("Node index not recognised.")
            
            if myStr[2] == 'x':
                component = 0
            elif myStr[2] == 'y':
                component = 1
            elif myStr[2] == 'z':
                component = 2
            else:
                raise IOError("Displacement component not recognised.")
            
            fp.write("{:<14,e}".format(PosDefor[index,component]))
            
        elif re.search(r'^M_[xyz]',myStr):
            if re.search(r'^M_[xyz]_[0-9]', myStr):
                index = int(myStr[4])
            elif re.search(r'root', myStr):
                index = 0
            elif re.search(r'tip', myStr):
                index = -1
            else:
                raise IOError("Node index not recognised.")
            
            if myStr[2] == 'x':
                component = 0
            elif myStr[2] == 'y':
                component = 1
            elif myStr[2] == 'z':
                component = 2
            else:
                raise IOError("Moment component not recognised.")
            
            locForces = localElasticForces(PosDefor,
                                           PsiDefor,
                                           PosIni,
                                           PsiIni,
                                           XBELEM,
                                           [index])
            
            fp.write("{:<14,e}".format(locForces[0,3+component]))
            
        elif re.search(r'^kappa_[xyz]',myStr):
            if re.search(r'^kappa_[xyz]_[0-9]', myStr):
                index = int(myStr[8])
            elif re.search(r'root', myStr):
                index = 0
            elif re.search(r'tip', myStr):
                index = -1
            else:
                raise IOError("Node index not recognised.")
            
            if myStr[6] == 'x':
                component = 0
            elif myStr[6] == 'y':
                component = 1
            elif myStr[6] == 'z':
                component = 2
            else:
                raise IOError("Moment component not recognised.")
            
            locStrain = localStrains(PosDefor,
                                     PsiDefor,
                                     PosIni,
                                     PsiIni,
                                     XBELEM,
                                     [index])
            
            fp.write("{:<14,e}".format(locStrain[0,3+component]))
        
        elif re.search(r'^u_.',myStr) and ctrlSurf != None:
            if re.search(r'^u_opt_[0-9]',myStr):
                index = int(myStr[6])
                fp.write("{:<14,e}".format(ctrlSurf.beta))
            else:
                raise IOError("u_opt output format not recognised")
            
        elif re.search(r'^du_.',myStr) and ctrlSurf != None:
            if re.search(r'^du_opt_[0-9]',myStr):
                index = int(myStr[7])
                fp.write("{:<14,e}".format(ctrlSurf.betaDot))
            else:
                raise IOError("du_opt output format not recognised")
            
        elif re.search(r'^contTime',myStr) and mpcCont != None:
            fp.write("{:<14,e}".format(mpcCont.contTime))
            
        else:
            ErrorMsg = "writeDict key not recognised. (" + myStr + ")"
            raise IOError(ErrorMsg)
    # END for myStr
    fp.write("\n")
    fp.flush()

if __name__ == '__main__':
    class Xbopts:
        Solution = ct.c_int(312)
    
    XBOPTS = Xbopts()
    writeDict = OrderedDict()
    writeDict['ya'] = 0
    OpenOutFile(writeDict, XBOPTS)