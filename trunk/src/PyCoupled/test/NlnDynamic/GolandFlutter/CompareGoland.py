'''
Created on 28 Feb 2013

@author: rjs10
'''

import numpy as np
import matplotlib.pyplot as plt

def TipHist():
    """@brief Tip history"""
    #        NumNodesElemArr = [3, 2]
    NumNodesElemArr = [2]
    
#        NumElemsArr = [12,24,36]
    NumElemsArr = [12,24,36]

    #M = [6, 12, 18]
    iMArr = [6]
    
#       U = [160.0, 161.0, 162.0, 163.0, 164.0, 165.0, 166.0, 167.0, 168.0, 169.0, \
#             170.0]
    UArr = [164.0]
    
    
    Rollup = False
    Tights = [False]
    KJMeth = True
    ImpStart = False
    
    "tuple of lines"
    lines = ()
    FFTs  = ()
    FFTfreqs = ()
    FFTlines = ()
    
    for NumNE in NumNodesElemArr:
        for NumElems in NumElemsArr:
            for iM in iMArr:
                for Tight in Tights:
                
                    "change to keep nodes equal"
                    if NumNE == 2:
                        NumElemsHere = NumElems
                    elif NumNE == 3:
                        NumElemsHere = int(0.5*NumElems)
                    else:
                        print("How many nodes per elem!!?")
                        return
                    
                    for U in UArr:
                        "read from file"
                        
                        Root = 'Convergence/'
                        
                        File = 'Goland' + \
                                '_NumNE' + str(NumNE) + \
                                '_NumEl' + str(NumElemsHere) + \
                                '_M' + str(iM) + '_U' + str(U)
                        
                        if Tight == True:
                            File += '_Tight'
                            
                        if Rollup == True:
                            File += '_Rollup'
                            
                        if KJMeth == True:
                            File += '_KJMeth'
                            
                        if ImpStart == True:
                            File += '_ImpStart'
                            
                        File += '_SOL312_dyn.dat'
                            
                        TipHist = np.loadtxt(Root + File \
                                             ,ndmin = 2, skiprows = 3)
                        
                        "tip z-displacement"
                        plt.figure(1)
                        line = plt.plot(TipHist[:,0],TipHist[:,3])
        #                if NumNE == 2:
        #                    plt.setp(line, marker = 'o', color = 'r')
        #                elif NumNE == 3:
        #                    plt.setp(line, marker = '+', color = 'k')
        #                
        #                if NumElems == 12:
        #                    plt.setp(line, ls = '-.')
        #                elif NumElems == 36:
        #                    plt.setp(line, ls = '-')
        
#                        if KJMeth == True:
#                            plt.setp(line, marker = 'o',ls = '-',color = 'k')
#                        elif KJMeth == False:
#                            plt.setp(line, marker = 's',ls = '-',color = 'r')
                        
                        lines = lines + (line,)
                        
                        "FFT"
                        dt = TipHist[1,0] - TipHist[0,0]
                        
                        FFTout = np.fft.fft(TipHist[:,3])
                        FFTfreq = np.fft.fftfreq(len(TipHist[:,3]),dt) * 2 * np.pi
                        
                        plt.figure(2)
                        FFTline = plt.plot(FFTfreq,abs(FFTout.imag))
                        
                        FFTs = FFTs + (FFTout,)
                        FFTfreqs = FFTfreqs + (FFTfreq,)
                        FFTlines = FFTlines + (FFTline,)
                    #End for U
                ##END for meth
            #END for M
        # END for NumElems
    #END for NumNodesElem
    
    plt.setp(lines, linewidth=2.0)
    plt.setp(lines, ms = 4.0)
    
    plt.figure(1)
    plt.xlabel('Physical time [s]')
    plt.ylabel('Vertical tip-deflection [m]')
#    plt.xlim(0,2.5)
#    plt.ylim(-2.0,2.0)
    plt.legend(("Loosely-coupled",\
                "Tightly-coupled"),\
                'lower left',shadow=True,prop={'size':12})
    
    plt.figure(2)
    plt.xlim(0.0,100.0)
    plt.xlabel('frequency [rads-1]')
    plt.ylabel('Frequency content [-] (abs(Im(FFT(z(t)))')
    
    plt.show()

if __name__ == '__main__':
    TipHist()
    
