import Main.SharPySettings as Settings
import PyBeam.Utils.DerivedTypes as DerivedTypes
import PyAero.UVLM.Utils.DerivedTypesAero as DerivedTypesAero
import PyCoupled.Utils.DerivedTypesAeroelastic as DerivedTypesAeroelastic
import ctypes as ct
import numpy as np
from PyBeam.Utils.XbeamLib import Skew
from collections import OrderedDict

# Beam options.
XBOPTS = DerivedTypes.Xbopts(FollowerForce = ct.c_bool(False),
                                 MaxIterations = ct.c_int(50),
                                 PrintInfo = ct.c_bool(True),
                                 NumLoadSteps = ct.c_int(25),
                                 Solution = ct.c_int(312),
                                 MinDelta = ct.c_double(1e-5),
                                 NewmarkDamp = ct.c_double(5e-3))
# beam inputs.
XBINPUT = DerivedTypes.Xbinput(3,14)
XBINPUT.BeamLength = 6.096
XBINPUT.BeamStiffness[0,0] = 1.0e+09
XBINPUT.BeamStiffness[1,1] = 1.0e+09
XBINPUT.BeamStiffness[2,2] = 1.0e+09
XBINPUT.BeamStiffness[3,3] = 0.99e+06
XBINPUT.BeamStiffness[4,4] = 9.77e+06
XBINPUT.BeamStiffness[5,5] = 1.0e+09
XBINPUT.BeamStiffness[:,:] = 1.0*XBINPUT.BeamStiffness[:,:]
XBINPUT.BeamMass[0,0] = 35.71
XBINPUT.BeamMass[1,1] = 35.71
XBINPUT.BeamMass[2,2] = 35.71
XBINPUT.BeamMass[3,3] = 8.64
XBINPUT.BeamMass[4,4] = 0.001
XBINPUT.BeamMass[5,5] = 0.001
# Off diagonal terms (in Theodorsen sectional coordinates)
ElasticAxis = -0.34
InertialAxis = -7.0/50.0
x_alpha = InertialAxis - ElasticAxis
# pitch-plunge coupling term (b-frame coordinates)
c = 1.8288
cgLoc = 0.5*c*np.array([0.0, x_alpha, 0.0])
cgSkew = Skew(cgLoc)
XBINPUT.BeamMass[:3,3:] = XBINPUT.BeamMass[0,0] * cgSkew
XBINPUT.BeamMass[3:,:3] = XBINPUT.BeamMass[:3,3:].T
XBINPUT.BeamMass[4,4] = 0.001 + XBINPUT.BeamMass[0,0]*pow(cgLoc[2],2.0)
XBINPUT.BeamMass[5,5] = 0.001 + XBINPUT.BeamMass[0,0]*pow(cgLoc[1],2.0)


# Get suggested panelling.
Umag = 140.0
M = 16
delTime = c/(Umag*M)

# Unsteady parameters.
XBINPUT.dt = delTime
XBINPUT.t0 = 0.0
XBINPUT.tfin = 0.5


# aero params.
WakeLength = 30.0*c
Mstar = int(WakeLength/(delTime*Umag))
# aero options.
N = XBINPUT.NumNodesTot - 1
VMOPTS = DerivedTypesAero.VMopts(M = M,
                                 N = N,
                                 ImageMethod = True,
                                 Mstar = Mstar,
                                 Steady = False,
                                 KJMeth = True,
                                 NewAIC = True,
                                 DelTime = delTime,
                                 NumCores = 4)
# Aero inputs.
# Control surface.
iMin = M - M/4
iMax = M
jMin = N - N/4
jMax = N
typeMotion = 'sin'
betaBar = 1.0*np.pi/180.0
omega = 15.0
ctrlSurf = DerivedTypesAero.ControlSurf(iMin,
                                        iMax,
                                        jMin,
                                        jMax,
                                        typeMotion,
                                        betaBar,
                                        omega)
# Inputs.
VMINPUT = DerivedTypesAero.VMinput(c = c,
                                   b = XBINPUT.BeamLength,
                                   U_mag = Umag,
                                   alpha = 0.0*np.pi/180.0,
                                   theta = 0.0,
                                   WakeLength = WakeLength,
                                   ctrlSurf = ctrlSurf)
# Unsteady aero inputs.
VelA_G   = np.array(([0.0,0.0,0.0]))
OmegaA_G = np.array(([0.0,0.0,0.0]))
VMUNST   = DerivedTypesAero.VMCoupledUnstInput(VMOPTS, VMINPUT, 0.0, 0.0,
                                               VelA_G, OmegaA_G)

# Aerolastic simulation results.
AELAOPTS = DerivedTypesAeroelastic.AeroelasticOps(ElasticAxis = ElasticAxis,
                                                  InertialAxis = InertialAxis,
                                                  AirDensity = 1.02,
                                                  Tight = False,
                                                  ImpStart = False)

# Live output options.
writeDict = OrderedDict()
writeDict['R_z (tip)'] = 0
writeDict['M_x (root)'] = 0
writeDict['M_y (root)'] = 0
writeDict['M_z (root)'] = 0