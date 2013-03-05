'''
Created on 25 Feb 2013

@author: rjs10
'''

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    
    "define Vels"
    Vels = np.array([160.0, 165.0, 170.0])
    
    "create tuple for holding line"
    lines = ()
    FFTs  = ()
    FFTfreqs = ()
    lines2 = ()
    FFTs2  = ()
    FFTfreqs2 = ()
    
    "create plot"
    plt.figure(1)
    
    "read data"
    for Vel in range(len(Vels)):
        TipHist = np.loadtxt('Goland_U' + str(Vels[Vel]) +'_SOL312_dyn.dat'\
                             ,ndmin = 2, skiprows = 3)
        
        TipHist2 = np.loadtxt('Goland_3noded_U' + str(Vels[Vel]) +'_SOL312_dyn.dat'\
                             ,ndmin = 2, skiprows = 3)
        
        line = plt.plot(TipHist[:,0],TipHist[:,3])
        line2 = plt.plot(TipHist2[:,0],TipHist2[:,3])
        
        dt = TipHist[1,0] - TipHist[0,0]
        dt2 = TipHist2[1,0] - TipHist2[0,0]
        
        FFTout = np.fft.fft(TipHist[:,3])
        FFTfreq = np.fft.fftfreq(len(TipHist[:,3]),dt) * 2 * np.pi
        FFTout2 = np.fft.fft(TipHist2[:,3])
        FFTfreq2 = np.fft.fftfreq(len(TipHist2[:,3]),dt) * 2 * np.pi
        
        lines = lines + (line,)
        FFTs = FFTs + (FFTout,)
        FFTfreqs = FFTfreqs + (FFTfreq,)
        lines2 = lines2 + (line2,)
        FFTs2 = FFTs2 + (FFTout2,)
        FFTfreqs2 = FFTfreqs2 + (FFTfreq2,)
    #END for Vel
    
    plt.setp(lines, linewidth=2.0)
    plt.setp(lines, marker = 'o', ms = 5.0)
    
    plt.xlabel('Physical time [s]')
    plt.ylabel('Vertical tip-deflection [m]')
    plt.xlim(0,1.0)
    plt.ylim(-0.005,0.005)
    plt.legend(("U = 160","U = 165","U = 170"),\
                'lower left',shadow=True,prop={'size':12})
    
    plt.figure(2)
    plt.plot(FFTfreqs[0], abs(FFTs[0].imag))
    plt.plot(FFTfreqs[1], abs(FFTs[1].imag))
    plt.plot(FFTfreqs[2], abs(FFTs[2].imag))
    plt.legend(("U = 160","U = 165","U = 170"),\
                'upper left',shadow=True,prop={'size':12})
    plt.xlim(0.0,100.0)
    plt.xlabel('frequency [rads-1]')
    plt.ylabel('Frequency content [-] (abs(Im(FFT(z(t)))')
    
    
    plt.show()
        