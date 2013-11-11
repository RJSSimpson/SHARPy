'''
@author: rjs10
'''
import numpy as np
import PyMPC.ssdiscrete

sizeFactor = 10

Ac = np.array([[-1.2822,       0.0,      0.98,     0.0],
                 [    0.0,       0.0,       1.0,     0.0],
                 [-5.4293,       0.0,   -1.8366,     0.0],
                 [ -128.2,     128.2,       0.0,     0.0]])


# 1 input
Bc = np.array([[ -0.3],
                 [  0.0],
                 [-17.0],
                 [  0.0]])

# 3 outputs
Cc = np.array([[    0.0,       1.0,       0.0,     0.0],
                 [    0.0,       0.0,       0.0,     1.0],
                 [ -128.2,     128.2,       0.0,     0.0]])

Dc = np.array([[  0.0],
                 [  0.0],
                 [  0.0]])

# Define state-space system.
dSys = PyMPC.ssdiscrete.StateSpace(Ac,Bc,Cc,Dc) # First create a continuous time.

dt = 0.5
dSys.cont2discrete(dt) # Convert to discrete.

# begin defining control using muaompc
N = 10
mu = 100

# TODO: Automatically append auxiliary states for slew rates.
AdSmall = np.zeros((5,5))
AdSmall[:AdSmall.shape[0]-1,:AdSmall.shape[1]-1] = np.array(dSys.A, copy=True)

BdSmall = np.array(dSys.B, copy=True)
BdSmall.resize((BdSmall.shape[0]+1, BdSmall.shape[1]))
BdSmall[dSys.nX:,:] = np.eye(dSys.nU)

# Make system large
Ad = np.zeros([AdSmall.shape[0]*sizeFactor,AdSmall.shape[0]*sizeFactor])
for i in range(sizeFactor):
    Ad[i*AdSmall.shape[0]:(i+1)*AdSmall.shape[0], i*AdSmall.shape[0]:(i+1)*AdSmall.shape[0]] = AdSmall
    
Bd = np.zeros([BdSmall.shape[0]*sizeFactor,1])
for j in range(sizeFactor):
    Bd[j*BdSmall.shape[0]:(j+1)*BdSmall.shape[0],:BdSmall.shape[1]] = BdSmall
    
# Cd and Dd not required for control

# Weighting matrices
Q = np.eye(Ad.shape[0])
Q[AdSmall.shape[0]-1::AdSmall.shape[0],AdSmall.shape[0]-1::AdSmall.shape[0]] = 0.0 #Auxiliary state for slew rates
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
KxSmall = np.array([[     0,      1,      0,      0,      0],
                    [-128.2,  128.2,      0,      0,      0],
                    [     0.,     0.,     0.,     0.,    -1.]])
Kx = np.zeros([KxSmall.shape[0],AdSmall.shape[1]*sizeFactor])
Kx[:KxSmall.shape[1],:AdSmall.shape[0]] = KxSmall
Ku = [[0],
      [0],
      [1]]

# Terminal constraints
f_lb = e_lb
f_ub = e_ub
F = Kx # Terminal constraint matrix
if __name__ == '__main__':
    pass