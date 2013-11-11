'''
Created on 24 Oct 2013

@author: rjs10
'''
import PyMPC.ssdiscrete as ssdiscrete
import numpy as np

# load reduced system from matlab

matPath = '/home/rjs10/git/SHARP/output/Goland/TorsionBending_M8N20_CS80/Q140_N50_Py'
disSys = ssdiscrete.StateSpace(matPath)

n = disSys.nX # Number of states.
m = disSys.nU # Number of inputs.
p = disSys.nY # Number of outputs.
dt = disSys.dt # model time step.

nAux = 1 # One auxiliary state for input rate constraint.

# Resize T to include zeros for auxiliary states.
T = np.zeros((n+nAux,disSys.T.shape[1]))
T[:n,:] = disSys.T.copy('C')

# Define for muaompc.

N = 100
mu = 100

Ad = np.zeros((n+nAux,n+nAux)) # MPC model state transition matrix.
Ad[:n,:n] = disSys.A.copy('C')

Bd = np.zeros((n+nAux,m)) # MPC model input matrix.
Bd[:n,:m] = disSys.B.copy('C')

if nAux > 0:
    Bd[n:,:] = np.eye(m)
    
# Weighting matrices
Q = 100.0*np.eye(Ad.shape[0])
Q[n:,n:] = 1.0 #Auxiliary state for slew rates
P = 10.0*Q[:,:] # Terminal constraint
R = 1.0*np.eye(Bd.shape[1])

# Input constraints
eui = 15.0*np.pi/180.0
u_lb = [[-eui]]
u_ub = [[eui]]

# State and output constraints
c = 0 # zero state and output constraints

# Mixed constraint matrices
Kx = np.zeros((c+nAux,n+nAux))
Ku = np.zeros((c+nAux,m))
e_lb = np.zeros((c+nAux,1))
e_ub = np.zeros((c+nAux,1))

if nAux > 0:
    Kx[-1] = -1.0
    Ku[-1] = 1.0
    inputRate = 30.0*np.pi/180.0*dt
    e_lb[-1,0] = -inputRate
    e_ub[-1,0] = inputRate
    
f_lb = e_lb
f_ub = e_ub
F = Kx # Terminal constraint matrix