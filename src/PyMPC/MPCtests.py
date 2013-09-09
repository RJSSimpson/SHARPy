'''
@author: Rob Simpson (rjs10@imperial.ac.uk)
http://ifatwww.et.uni-magdeburg.de/syst/research/muAO-MPC/doc/html/tutor.advanced.html
'''

from muaompc.ltidt import setup_mpc_problem
import numpy as np
import sys, os
import cvxopt

if __name__ == '__main__':
    
    sys.path.append(os.getcwd() + '/systems')
    
    # set-up MPC problem
    mpc = setup_mpc_problem('sys_aircraftLargeSys')
    mpc.generate_c_files()
    
    # Configure optimizer
    ctl = mpc.ctl
    ctl.conf.in_iter = 24 # Internal iterations in augmented Lagrangian method
    ctl.conf.ex_iter  =2  # External iteration in augmented Lagrangian method
    ctl.conf.warmstart = True
    
    # Initialize state
    x = np.zeros([mpc.size.states, 1])
    x[3] = 0.0
    
    # Set point
    setPoint = np.zeros([mpc.size.states,1])
    # setPoint[3] = 400.0
    setPoint[3::5] = 400.0
    
    for i in range(mpc.N):
        ctl.x_ref[mpc.size.states*i:mpc.size.states*(i+1)] = setPoint
        
    # form QP data at Initial condition
    ctl.form_qp(x)
    

        
    # Number of timesteps
    s = 40
    for i in range(s):
        ctl.solve_problem(x)
        print("step:", i, "\t u[0] = ", ctl.u_opt[0], "\t x = ", x[:4].T)
        x = ctl.sys.Ad.dot(x) + ctl.sys.Bd.dot(ctl.u_opt[:mpc.size.inputs])
    print("SUCCESS!")