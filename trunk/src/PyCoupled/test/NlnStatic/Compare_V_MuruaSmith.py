'''
Created on 15 Feb 2013

@author: rjs10
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#matplotlib with latex
from matplotlib import rc
rc('font',**{'family':'serif','sans-serif':['Times New Roman']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

mpl.rcParams.update({'font.size': 20})

if __name__ == '__main__':
    
#    Sharp2 = np.loadtxt('Murua2012_Smith2001_SHARP2.dat',ndmin = 2)
#    Sharp4 = np.loadtxt('Murua2012_Smith2001_SHARP4.dat',ndmin = 2)
    Smith2 = np.loadtxt('Murua2012_Smith2001_SMITH2.dat',ndmin = 2)
    Smith4 = np.loadtxt('Murua2012_Smith2001_SMITH4.dat',ndmin = 2)
    
    plt.figure(1)
    plt.xlim((0.0,16.0))
#    plt.ylim((0.0, 6.0))
    plt.axis('equal')
#    a = plt.plot(Sharp2[:,0],Sharp2[:,1],'k-')
    a = plt.plot(Smith2[:,0],Smith2[:,1],'k--')
#    c = plt.plot(Sharp4[:,0],Sharp4[:,1],'r-')
    b = plt.plot(Smith4[:,0],Smith4[:,1],'r--')
    
    SHARPy2 = np.loadtxt('CantHALE2_SOL112_def.dat',skiprows=2,ndmin = 2)
    SHARPy4 = np.loadtxt('CantHALE4_SOL112_def.dat',skiprows=2,ndmin = 2)
    
    c = plt.plot(SHARPy2[:,2],SHARPy2[:,4],'ko-')
    d = plt.plot(SHARPy4[:,2],SHARPy4[:,4],'ro-')
    lines = a,b,c,d
    plt.setp(lines,linewidth=2,markersize=6)
    plt.xlabel(r'Spanwise ordinate [$m$]')
    plt.ylabel('Vertical deflection [$m$]')
    plt.legend((r'Smith et al.(2001), $\alpha$ = $2$ \textit{deg}',\
                r'Smith et al.(2001), $\alpha$ = $4$ \textit{deg}',\
                r'Current, $\alpha$ = $2$ \textit{deg}',\
                r'Current, $\alpha$ = $4$ \textit{deg}'),\
               'upper left',shadow=False,prop={'size':20})
    plt.grid(True)
    plt.show()
    
    
    "plot convergence"
    KJ = np.loadtxt('ConvergenceKJ.dat',ndmin = 2)
    KP = np.loadtxt('ConvergenceKP.dat',ndmin = 2)
    KJ_noTip = np.loadtxt('ConvergenceKJ_KJnoTip.dat',ndmin = 2)
    
    plt.figure(2)
    g = plt.plot([1,  5,  10,  15,  20,  30],KJ[:,4])
    h = plt.plot([1,  5,  10,  15,  20,  30],KP[:,4])
    plt.xlabel('3-Noded beam elements')
    plt.ylabel('Tip displacement [m]')
    plt.legend(("Joukowski method","Katz and Plotkin method"),'lower right',shadow=True,prop={'size':12})
    plt.show()
    
    plt.figure(3)
    g = plt.plot([1,  2,  5,  10,  15],KJ[4,:])
    h = plt.plot([1,  2,  5,  10,  15],KP[4,:])
    plt.xlabel('Chordwise panels, $M$')
    plt.ylabel('Tip displacement [m]')
    plt.legend(("Joukowski method","Katz and Plotkin method"),'lower right',shadow=True,prop={'size':12})
    plt.show()
    
#    plt.figure(4)
#    g = plt.plot([1,  2,  5,  10,  15],KJ[4,:])
#    h = plt.plot([1,  2,  5,  10,  15],KJ_noTip[4,:])
#    plt.xlabel('Chordwise Panels')
#    plt.ylabel('Tip displacement [m]')
#    plt.legend(("Joukowski Method","Without tip contribution"),'lower right',shadow=True,prop={'size':12})
#    plt.show()
    
    