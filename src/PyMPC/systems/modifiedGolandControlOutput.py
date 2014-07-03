'''
Created on 21 Nov 2013

@author: rjs10
'''
import PyMPC.ssdiscrete as ssdiscrete
import numpy as np
import ctypes as ct

# load reduced system from matlab

matPath = '/home/rjs10/git/SHARP/output/ModifiedGoland/TorsionBending_M8N20_CS80/Q28_N8_Py'
disSys = ssdiscrete.StateSpace(matPath)

n = disSys.nX # Number of states.
m = disSys.nU # Number of inputs.
p = disSys.nY # Number of outputs.
dt = disSys.dt # model time step.

nAux = 0 # One auxiliary state for input rate constraint.

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

Cd = np.zeros((p+nAux,n+nAux)) # used for output weighting only.
Cd[:p,:n] = disSys.C.copy('C')
if nAux > 0:
    Cd[-nAux:-nAux] = np.eye(nAux) # keep input rates as outputs

Dd = np.zeros((p+nAux,m)) # no feed-through, used for output weighting only.
Dd[:p,:] = disSys.D.copy('C')

if nAux > 0:
    Bd[n:,:] = np.eye(m)
    
# Define desired output weighting
Q_out = np.eye(p+nAux)
Q_out[0,0] = 0.0 # bending only
if nAux > 0:
    Q_out[-nAux:,-nAux:] = np.zeros((nAux,nAux))
R_out = np.eye(m)
    
# Weighting matrices based on input weighting
Q = np.dot(Cd.T,np.dot(Q_out,Cd))
R = np.dot(Dd.T,np.dot(Q_out,Dd)) + R_out
Q = Q/np.linalg.norm(Q) # normalize to give scale as with state weights
R = R/np.linalg.norm(R)
P = 1.0*Q[:,:] # Terminal constraint

# Input constraints
eui = 10.0*np.pi/180
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
#     Kx[-1,-1] = -1.0
#     Ku[-1,-1] = 1.0
#     inputRate = 300.0*np.pi/180.0*dt
#     e_lb[-1,0] = -inputRate
#     e_ub[-1,0] = inputRate
    pass
    
f_lb = e_lb
f_ub = e_ub
F = Kx #np.zeros_like(Kx,ct.c_double,'C') # Terminal constraint matrix