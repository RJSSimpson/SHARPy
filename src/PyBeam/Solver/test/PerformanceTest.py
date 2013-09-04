'''PyBeam.Solver.test.PerformanceTest
@brief      Timing of PyBeam Solvers.
@author     Rob Simpson
@contact    r.simpson11@imperial.ac.uk
@version    0.0
@date       23/11/2012
@pre        None
@warning    None
'''
import NonlinearStatic
import NonlinearDynamic
import DerivedTypes
from Misc import Timer
import numpy as np
import matplotlib.pyplot as plt

def TimeStaticSolvers():
        
        """Set up Xbopts for nonlinear static analysis defined in 
        NonlinearStatic/testcases.pdf case 1.1."""
        XBOPTS = DerivedTypes.Xbopts()
        XBOPTS.Solution.value = 112 
        XBOPTS.NumLoadSteps.value = 10
        XBOPTS.MinDelta.value = 1e-05
        XBOPTS.FollowerForce.value = False
        
        
        "Initialise timing result array"
        Times = np.zeros((10,4))
              
        Counter = 0
        for NumElems in range(10,110,10):
            
            """Set up Xbinput for nonlinear static analysis defined in 
            NonlinearStatic/testcases.pdf case 1.1."""
            XBINPUT = DerivedTypes.Xbinput(2,NumElems)
            XBINPUT.BeamLength = 5.0
            XBINPUT.BeamStiffness[0,0] = 4.8e+08
            XBINPUT.BeamStiffness[1,1] = 3.231e+08
            XBINPUT.BeamStiffness[2,2] = 3.231e+08
            XBINPUT.BeamStiffness[3,3] = 1.0e+06
            XBINPUT.BeamStiffness[4,4] = 9.346e+06
            XBINPUT.BeamStiffness[5,5] = 9.346e+06
            XBINPUT.BeamMass[0,0] = 100
            XBINPUT.BeamMass[1,1] = 100
            XBINPUT.BeamMass[2,2] = 100
            XBINPUT.BeamMass[3,3] = 10
            XBINPUT.BeamMass[4,4] = 0.0 #Neglect the cross-section bending inertia
            XBINPUT.BeamMass[5,5] = 0.0 #Neglect the cross-section bending inertia
            XBINPUT.ForceStatic[-1,2] = 6e+05
            
            print('Number of Elements = %s' %(NumElems))
            
            with Timer() as t_F90:
                NonlinearStatic.Solve_F90(XBINPUT, XBOPTS)
            
            print('Request took %.03f sec.' % t_F90.interval)
            
            
            with Timer() as t_Steps:
                NonlinearStatic.Solve_F90_steps(XBINPUT, XBOPTS)
                
            print('Request took %.03f sec.' % t_Steps.interval)
            
                
            with Timer() as t_Py:
                NonlinearStatic.Solve_Py(XBINPUT, XBOPTS)
                
            print('Request took %.03f sec.' % t_Py.interval)
            
            Times[Counter,:] = [NumElems, t_F90.interval,\
                                 t_Steps.interval, t_Py.interval]
            
            Counter += 1
        
        print(Times)
        plt.figure(1)
        plt.plot(Times[:,0],Times[:,1],'ko-')
        plt.plot(Times[:,0],Times[:,2],'rs-')
        plt.plot(Times[:,0],Times[:,3],'g^-')
        plt.xlabel('NumElems (2-noded)')
        plt.ylabel('Time, s.')
        plt.legend(('F90 Solver from Lib',\
                   'F90 w/ steps in Python',\
                   'Python Solver'),\
                   'upper left')
        plt.show()

    
def TimeDynamicSolvers():
    "Set-up problem"
        
    """Set up Xbopts for nonlinear static/dynamic analysis defined in 
    NonlinearStatic/testcases.pdf case 1.2."""
    XBOPTS = DerivedTypes.Xbopts()
    XBOPTS.Solution.value = 312 
    XBOPTS.NumLoadSteps.value = 10
    XBOPTS.MinDelta.value = 1e-05
    XBOPTS.FollowerForce.value = False
    XBOPTS.NumGauss.value = 1
    XBOPTS.PrintInfo.value = False
    
    "Initialise timing result array"
    Times = np.zeros((4,4))
          
    Counter = 0
    for NumElems in range(10,50,10):
        
        """Set up Xbinput for nonlinear static/dynamic analysis defined in 
        NonlinearStatic/testcases.pdf case 1.2."""
        XBINPUT = DerivedTypes.Xbinput(2,NumElems)
        XBINPUT.BeamLength = 5.0
        XBINPUT.BeamStiffness[0,0] = 4.8e+08
        XBINPUT.BeamStiffness[1,1] = 3.231e+08
        XBINPUT.BeamStiffness[2,2] = 3.231e+08
        XBINPUT.BeamStiffness[3,3] = 1.0e+06
        XBINPUT.BeamStiffness[4,4] = 9.346e+06
        XBINPUT.BeamStiffness[5,5] = 9.346e+06
        XBINPUT.BeamMass[0,0] = 100
        XBINPUT.BeamMass[1,1] = 100
        XBINPUT.BeamMass[2,2] = 100
        XBINPUT.BeamMass[3,3] = 10
        XBINPUT.BeamMass[4,4] = 0.0001
        XBINPUT.BeamMass[5,5] = 0.0001
        "Dynamic parameters"
        XBINPUT.t0 = 0.0
        XBINPUT.tfin = 2.0
        XBINPUT.dt = 0.001
        XBINPUT.Omega = 20.0
        XBINPUT.ForceDyn[-1,2] = 60e+03
        XBINPUT.ForcingType = 'RampSin'
        XBINPUT.RampTime = 1.0
        
        print('Number of Elements = %s' %(NumElems))
        
        with Timer() as t_F90:
            NonlinearDynamic.Solve_F90(XBINPUT, XBOPTS)
        
        print('Request took %.03f sec.' % t_F90.interval)
        
        
        with Timer() as t_Steps:
            NonlinearDynamic.Solve_F90_steps(XBINPUT, XBOPTS)
            
        print('Request took %.03f sec.' % t_Steps.interval)
        
            
        with Timer() as t_Py:
            NonlinearDynamic.Solve_Py(XBINPUT, XBOPTS)
            
        print('Request took %.03f sec.' % t_Py.interval)
        
        Times[Counter,:] = [NumElems, t_F90.interval,\
                             t_Steps.interval, t_Py.interval]
        
        Counter += 1
    
    print(Times)
    plt.figure(1)
    plt.plot(Times[:,0],Times[:,1],'ko-')
    plt.plot(Times[:,0],Times[:,2],'rs-')
    plt.plot(Times[:,0],Times[:,3],'g^-')
    plt.xlabel('NumElems (2-noded)')
    plt.ylabel('Time, s.')
    plt.legend(('F90 Solver from Lib',\
               'F90 w/ steps in Python',\
               'Python Solver'),\
               'upper left')
    plt.show()
    

if __name__ == '__main__':
    #TimeStaticSolvers()
    TimeDynamicSolvers()
    
    