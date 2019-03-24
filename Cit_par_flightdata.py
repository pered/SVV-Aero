# Citation 550 - Linear simulation



# Stationary flight condition
import scipy as np

#hp0    = 8000.      	     # pressure altitude in the stationary flight condition [m]
#V0     = 59.9          # true airspeed in the stationary flight condition [m/sec]
#alpha0 = 3. * np.pi / 180            # angle of attack in the stationary flight condition [rad]
#th0    = 2. * np.pi / 180            # pitch angle in the stationary flight condition [rad]

# Aircraft mass
#m      = 4547.8            # mass [kg]

# aerodynamic properties
e      = 0.9336            # Oswald factor [ ]
CD0    = 0.0215            # Zero lift drag coefficient [ ]
CLa    = 0.0771 * 180/np.pi       # Slope of CL-alpha curve [ ]

# Longitudinal stability
Cma    = -0.60211        # longitudinal stabilty [ ]
Cmde   = -1.4124           # elevator effectiveness [ ]

# Aircraft geometry

S      = 24.2          # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	          # [ ]
lh     = 5.5   # tail length [m]
c      = 2.022	          # mean aerodynamic cord [m]
lh_c   = lh / c	          # [ ]
b      = 15.911	          # wing span [m]
bh     = 5.791	          # stabilser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
Vh_V   = 1	          # [ ]
ih     = -2 * np.pi / 180   # stabiliser angle of incidence [rad]
xcg    = 0.3 * c

# Constant values concerning atmosphere and gravity

rho0   = 1.2250          # air density at sea level [kg/m^3] 
lambdaa = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)

# air density [kg/m^3]  
#rho    = rho0 * np.power( ((1+(lambdaa * hp0 / Temp0))), (-((g / (lambdaa*R)) + 1)))   
#W      = m * g            # [N]       (aircraft weight)

# Constant values concerning aircraft inertia

muc    = 102.7
#mub    = m / (rho * S * b)
KX2    = 0.019
KZ2    = 0.042
KXZ    = 0.002
KY2    = 0.980

# Aerodynamic constants

Cmac   = 0.                      # Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * np.pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Lift and drag coefficient

#CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]
#CD = CD0 + (CLa * alpha0) ** 2 / (np.pi * A * e) # Drag coefficient [ ]

# Stabiblity derivatives

#CX0    = 0.
CXu    = -0.02792
CXa    = -0.47966
CXadot = +0.0833
CXq    = -0.2817
CXde   = -0.0373
CXdt   = 0.0

#CZ0    = -1.136
CZu    = -0.3762
CZa    = -5.7434
CZadot = -0.0035
CZq    = -5.6629
CZde   = -0.6961
CZdt   = 0.0

Cmu    = +0.0699
Cmadot = -0.178
Cmq    = -8.7941
Cmde   = -1.1642
Cmdt   = 0.0

CYb    = -0.7500
CYbdot =  0.0     
CYp    = -0.0304
CYr    = +0.8495
CYda   = -0.0400
CYdr   = +0.2300

Clb    = -0.10260
Clp    = -0.71085
Clr    = +0.23760
Clda   = -0.23088
Cldr   = +0.03440

Cnb    =  +0.1348
Cnbdot =   0.0     
Cnp    =  -0.0602
Cnr    =  -0.2061
Cnda   =  -0.0120
Cndr   =  -0.0939
