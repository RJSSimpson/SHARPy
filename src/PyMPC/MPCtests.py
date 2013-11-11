'''
@author: Rob Simpson (rjs10@imperial.ac.uk)
http://ifatwww.et.uni-magdeburg.de/syst/research/muAO-MPC/doc/html/tutor.advanced.html
'''

from muaompc.ltidt import setup_mpc_problem
import numpy as np
import sys, os
import cvxopt

if __name__ == '__main__':
    
    sys.path.append(os.getcwd() + '/systems/')
    
    # set-up MPC problem
    mpc = setup_mpc_problem('sys_aircraft')
    mpc.generate_c_files()
    
    # Configure optimizer
    mpc.ctl.conf.in_iter = 24 # Internal iterations in augmented Lagrangian method
    mpc.ctl.conf.ex_iter  = 2  # External iteration in augmented Lagrangian method
    mpc.ctl.conf.warmstart = True
    
    # Initialize state
    x = np.zeros([mpc.size.states, 1])
    x[3] = 0.0
    
    # Set point
    setPoint = np.zeros([mpc.size.states,1])
    setPoint[3] = 400.0
    
    print("set point = ",setPoint)
    
    for i in range(mpc.N):
        mpc.ctl.x_ref[mpc.size.states*i:mpc.size.states*(i+1)] = setPoint
        
    print("x_ref = ",mpc.ctl.x_ref)
    print("shape x_ref",np.shape(mpc.ctl.x_ref))
        
    # form QP data at Initial condition
    mpc.ctl.form_qp(x)
       
    # Number of timesteps
    s = 100
    for i in range(s):
        
        # calculate disturbance
        w = [[0.0],
             [0.0],
             [0.0],
             [1.0],
             [0.0]] #np.random.randn()
        
        if i == 50:
            setPoint[3] = 800 
            print("Set point change = ",setPoint)
        for j in range(mpc.N):
            mpc.ctl.x_ref[mpc.size.states*j:mpc.size.states*(j+1)] = (setPoint - w)
        # solve MPC problem.
        
        mpc.ctl.solve_problem(x)
        
        print("step:", i, "\t u[0] = ", mpc.ctl.u_opt[0], "\t x = ", x[:4].T)
        x = mpc.ctl.sys.Ad.dot(x) + mpc.ctl.sys.Bd.dot(mpc.ctl.u_opt[:mpc.size.inputs]) + w
    print("FINISHED!")