'''
@author: rjs10
'''
import numpy as np
import PyMPC.ssdiscrete as ssdiscrete

A = np.array([[-1.2822,       0.0,      0.98,     0.0],
                 [    0.0,       0.0,       1.0,     0.0],
                 [-5.4293,       0.0,   -1.8366,     0.0],
                 [ -128.2,     128.2,       0.0,     0.0]])

# 1 input
B = np.array([[ -0.3],
                 [  0.0],
                 [-17.0],
                 [  0.0]])

# 3 outputs
C = np.array([[    0.0,       1.0,       0.0,     0.0],
                 [    0.0,       0.0,       0.0,     1.0],
                 [ -128.2,     128.2,       0.0,     0.0]])

D = np.array([[  0.0],
                 [  0.0],
                 [  0.0]])

# Define state-space system.
dSys = ssdiscrete.StateSpace(A,B,C,D) # First create a continuous time.

dt = 0.5
dSys.cont2discrete(dt) # Convert to discrete.

# begin defining control using muaompc
N = 10
mu = 100

# TODO: Automatically append auxiliary states for slew rates.
Ad = np.zeros((5,5))
Ad[:Ad.shape[0]-1,:Ad.shape[1]-1] = np.array(dSys.A, copy=True)

Bd = np.array(dSys.B, copy=True)
Bd.resize((Bd.shape[0]+1, Bd.shape[1]))
Bd[dSys.nX:,:] = np.eye(dSys.nU)

# Weighting matrices
Q = 10.0*np.eye(Ad.shape[0])
Q[4,4] = 0.0 #Auxiliary state for slew rates
P = Q # Terminal constraint
R = np.eye(Bd.shape[1])

# Input constraints
eui = 15.0*np.pi/180.0
u_lb = [[-eui]]
u_ub = [[eui]]

# Mixed constraints
ex2 = 20.0*np.pi/180.0 # Pitch angle.
ex5 = 30*np.pi/180.0*dt # Maximum change in control input each time step.
ey3 = 30.0 # Altitude rate limit, ms-1.
e_lb = [[-ex2], [-ey3], [-ex5]]
e_ub = [[ex2], [ey3], [ex5]]

# Mixed constraint matrices
Kx = [[     0,      1,      0,      0,      0],
      [-128.2,  128.2,      0,      0,      0],
      [     0.,     0.,     0.,     0.,    -1.]]
Ku = [[0],
      [0],
      [1]]

# Terminal constraints
f_lb = e_lb
f_ub = e_ub
F = Kx # Terminal constraint matrix