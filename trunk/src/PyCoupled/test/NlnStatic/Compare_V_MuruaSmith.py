'''
Created on 15 Feb 2013

@author: rjs10
'''

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    
    "plot versus smith"
    dt = np.dtype([('Span',float),\
                   ('Vert',float)])
    
    Sharp2 = np.loadtxt('Murua2012_Smith2001_SHARP2.dat',ndmin = 2)
    Sharp4 = np.loadtxt('Murua2012_Smith2001_SHARP4.dat',ndmin = 2)
    Smith2 = np.loadtxt('Murua2012_Smith2001_SMITH2.dat',ndmin = 2)
    Smith4 = np.loadtxt('Murua2012_Smith2001_SMITH4.dat',ndmin = 2)
    
    plt.figure(1)
    plt.axis('equal')
    plt.xlim((0.0,16.0))
    plt.ylim((0.0, 6.0))
    a = plt.plot(Sharp2[:,0],Sharp2[:,1],'k-')
    b = plt.plot(Smith2[:,0],Smith2[:,1],'k--')
    c = plt.plot(Sharp4[:,0],Sharp4[:,1],'r-')
    d = plt.plot(Smith4[:,0],Smith4[:,1],'r--')
    
    SHARPy2 = np.loadtxt('CantHALE2_SOL112_def.dat',skiprows=2,ndmin = 2)
    SHARPy4 = np.loadtxt('CantHALE4_SOL112_def.dat',skiprows=2,ndmin = 2)
    
    e = plt.plot(SHARPy2[:,2],SHARPy2[:,4],'ko')
    f = plt.plot(SHARPy4[:,2],SHARPy4[:,4],'ro')
    lines = a,b,c,d,e,f
    plt.setp(lines,linewidth=2,markersize=6)
    plt.xlabel('Spanwise ordinate [m]')
    plt.ylabel('Vertical deflection [m]')
    plt.legend(("SHARP - alpha = 2 deg","Smith et al. 2001 - alpha = 2 deg",\
                "SHARP - alpha = 4 deg","Smith et al. 2001 - alpha = 4 deg",\
                "SHARPy - alpha = 2 deg", "SHARPy - alpha = 4 deg"),\
               'upper left',shadow=True,prop={'size':12})
    plt.show()
    
    
    "plot convergence"
    KJ = np.loadtxt('ConvergenceKJ.dat',ndmin = 2)
    KP = np.loadtxt('ConvergenceKP.dat',ndmin = 2)
    KJ_noTip = np.loadtxt('ConvergenceKJ_KJnoTip.dat',ndmin = 2)
    
    plt.figure(2)
    g = plt.plot([1,  5,  10,  15,  20,  30],KJ[:,4])
    h = plt.plot([1,  5,  10,  15,  20,  30],KP[:,4])
    plt.xlabel('3-Noded Beam Elements')
    plt.ylabel('Tip displacement [m]')
    plt.legend(("Joukowski Method","Katz and Plotkin Method"),'lower right',shadow=True,prop={'size':12})
    plt.show()
    
    plt.figure(3)
    g = plt.plot([1,  2,  5,  10,  15],KJ[4,:])
    h = plt.plot([1,  2,  5,  10,  15],KP[4,:])
    plt.xlabel('Chordwise Panels')
    plt.ylabel('Tip displacement [m]')
    plt.legend(("Joukowski Method","Katz and Plotkin Method"),'lower right',shadow=True,prop={'size':12})
    plt.show()
    
#    plt.figure(4)
#    g = plt.plot([1,  2,  5,  10,  15],KJ[4,:])
#    h = plt.plot([1,  2,  5,  10,  15],KJ_noTip[4,:])
#    plt.xlabel('Chordwise Panels')
#    plt.ylabel('Tip displacement [m]')
#    plt.legend(("Joukowski Method","Without tip contribution"),'lower right',shadow=True,prop={'size':12})
#    plt.show()
    
    