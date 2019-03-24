# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:27:49 2019

@author: M.Osowski
"""

from Cit_par_flightdata import *
import scipy.io
import scipy as np
import matplotlib.pyplot as plt
from Cit_class_experimental import Symmetrical_Model_Numerical#, Asymmetric_Model_numerical
from Flight_data_reader import plot_data


 
    
Flightdata = scipy.io.loadmat('FTISxprt-20190319_123022.mat')

TimeData = Flightdata['flightdata']['time'][0][0][0][0][0].transpose()

VaneAOA = Flightdata['flightdata']['vane_AOA'][0][0][0][0][0]
Pitch = Flightdata['flightdata']['Ahrs1_Pitch'][0][0][0][0][0]
Psi = Flightdata['flightdata']['Ahrs1_Roll'][0][0][0][0][0]

d_e = Flightdata['flightdata']['delta_e'][0][0][0][0][0]
d_et = Flightdata['flightdata']['elevator_dte'][0][0][0][0][0]

VTAS = Flightdata['flightdata']['Dadc1_tas'][0][0][0][0][0] * 0.514444 #converted to m/s
hp = Flightdata['flightdata']['Dadc1_alt'][0][0][0][0][0] * 0.3048 #converted to m

Fuel_used_l = Flightdata['flightdata']['lh_engine_FU'][0][0][0][0][0] * 0.453592
Fuel_used_r = Flightdata['flightdata']['rh_engine_FU'][0][0][0][0][0] * 0.453592
Fuel_used = Fuel_used_l + Fuel_used_r
m0 = 14730. * 0.453592

def compare_s(tstart, tend):
    
    #** set parameters **
    a0 = VaneAOA[(TimeData >= tstart) & (TimeData <= tend)][0]
    theta0 = Pitch[(TimeData >= tstart) & (TimeData <= tend)][0]
    V_0 = VTAS[(TimeData >= tstart) & (TimeData <= tend)][0]
    hp0 = hp[(TimeData >= tstart) & (TimeData <= tend)][0] 
    m = m0 - Fuel_used[(TimeData >= tstart) & (TimeData <= tend)][0]
    rho = rho0 * np.power( ((1+(lambdaa * hp0 / Temp0))), (-((g / (lambdaa*R)) + 1))) 
    W = m * g
    mub = m / (rho * S * b)
    CL = 2 * W / (rho * V_0 ** 2 * S)
    CD = CD0 + (CLa * a0) ** 2 / (np.pi * A * e)
    CX0 = CL * np.sin(theta0)
    CZ0 = -CL * np.cos(theta0)
    print(CX0,CZ0)

    
    timedomain = TimeData[(TimeData >= tstart) & (TimeData <= tend)]
    de = d_e[(TimeData >= tstart) & (TimeData <= tend)] 
    print(max(de))
    det = d_et[(TimeData >= tstart) & (TimeData <= tend)]
    de_in = np.vstack((de, det)) * np.pi/180.
    
    timeinterval = timedomain[-1] - timedomain[0]
    
    model = Symmetrical_Model_Numerical(V_0, c, CXu, CXa, CZ0, CZu, CZa, CZadot, muc, CZq, Cmu, Cmadot, KY2, Cma, CX0, Cmq, CXde, CXdt, CZde, CZdt, Cmde, Cmdt)
    plt.figure()
#    model.pulse()
    model.forced(de_in, timeinterval, 0.1, 0, 0, 0, 0.)
#    plt.figure()
#    model.eigs()
#    plt.plot(timedomain - tstart, de)
    
    plot_data(tstart, tend, [Pitch-theta0, VTAS])

def compare_a(tstart, tend):

    #** set parameters **
    b0 = 0.
    psi0 = Psi[(TimeData >= tstart) & (TimeData <= tend)][0]
    
    
    




compare_s(2780, 3050)


