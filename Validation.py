# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:27:49 2019

@author: M.Osowski
"""

from Cit_par_1 import *
import scipy.io
import scipy as np
import matplotlib.pyplot as plt
from Cit_class_experimental import Symmetrical_Model_Numerical, Asymmetric_Model_numerical
from Flight_data_reader import plot_data

plt.close('all')
    
Flightdata = scipy.io.loadmat('FTISxprt-20190319_123022.mat')

TimeData = Flightdata['flightdata']['time'][0][0][0][0][0].transpose()

VaneAOA = Flightdata['flightdata']['vane_AOA'][0][0][0][0][0]
Pitch = Flightdata['flightdata']['Ahrs1_Pitch'][0][0][0][0][0]
Psi = Flightdata['flightdata']['Ahrs1_Roll'][0][0][0][0][0]
q = Flightdata['flightdata']['Ahrs1_bPitchRate'][0][0][0][0][0]
p = Flightdata['flightdata']['Ahrs1_bRollRate'][0][0][0][0][0]
r = Flightdata['flightdata']['Ahrs1_bYawRate'][0][0][0][0][0]

d_e = Flightdata['flightdata']['delta_e'][0][0][0][0][0]
d_et = Flightdata['flightdata']['elevator_dte'][0][0][0][0][0]
d_a = Flightdata['flightdata']['delta_a'][0][0][0][0][0]
d_r = Flightdata['flightdata']['delta_r'][0][0][0][0][0]

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
    q0 = q[(TimeData >= tstart) & (TimeData <= tend)][0]
    
    hp0 = hp[(TimeData >= tstart) & (TimeData <= tend)][0] 
    m = m0 - Fuel_used[(TimeData >= tstart) & (TimeData <= tend)][0]
    rho = rho0 * np.power( ((1+(Lambda * hp0 / Temp0))), (-((g / (Lambda*R)) + 1))) 
    
    W = m * g
    muc = m / (rho * S * c)
    print(V_0)
    CL = 2 * W / (rho * V_0 ** 2 * S)
    CD = CD0 + (CLa * a0) ** 2 / (np.pi * A * e)
    CX0 = CL * np.sin(theta0*np.pi/180.)
    CZ0 = -CL * np.cos(theta0*np.pi/180.)
    print(CX0, CZ0)
    q_in = q0*c/V_0
    
    timedomain = TimeData[(TimeData >= tstart) & (TimeData <= tend)]
    de = d_e[(TimeData >= tstart) & (TimeData <= tend)] - d_e[(TimeData >= tstart) & (TimeData <= tend)][0]
    det = d_et[(TimeData >= tstart) & (TimeData <= tend)]
    de_in = np.vstack((de, det)) * np.pi/180.
    
    timeinterval = timedomain[-1] - timedomain[0]

    model = Symmetrical_Model_Numerical(V_0, c, CXu, CXa, CZ0, CZu, CZa, CZadot, muc, CZq, Cmu, Cmadot, KY2, Cma, CX0, Cmq, CXde, CXdt, CZde, CZdt, Cmde, Cmdt)
    
    tt_model, out_model = model.forced(de_in, timeinterval, 0.1, 0, 0, 0, 0 , hp0)
    
    
    #**plotting**
    
    plt.figure()
    plt.figlegend('Model response')
    
    plt.subplot(231) #x Velocity
    plt.title('x velocity')
    plt.plot(timedomain, out_model[0], label = 'Model response')
    plt.plot(timedomain, VTAS[(TimeData >= tstart) & (TimeData <= tend)])
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity [m/s]')
    
    plt.subplot(232) #AOA
    plt.title('Angle of attack')
    plt.plot(timedomain, out_model[1])
    plt.plot(timedomain, VaneAOA[(TimeData >= tstart) & (TimeData <= tend)] - a0)
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [°]')
    
    plt.subplot(233) #Pitch
    plt.title('Pitch attitude')
    plt.plot(timedomain, out_model[2])
    plt.plot(timedomain, Pitch[(TimeData >= tstart) & (TimeData <= tend)] - theta0)
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [°]')
    
    plt.subplot(234) #q
    plt.title('Pitch rate')
    plt.plot(timedomain, out_model[3])
    plt.plot(timedomain, q[(TimeData >= tstart) & (TimeData <= tend)]- q0)
    plt.xlabel('Time [s]')
    plt.ylabel('Angular velocity [°/s]')
    
    plt.subplot(235) #Altitude
    plt.title('Altitude')
    plt.plot(timedomain, out_model[4])
    plt.plot(timedomain, hp[(TimeData >= tstart) & (TimeData <= tend)])
    plt.xlabel('Time [s]')
    plt.ylabel('Altitude [m]')
    
    plt.subplot(236) #Elevator
    plt.title('Measured elevator deflection')
    plt.plot(timedomain, de, 'darkorange')
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [°]')
    
    plt.figlegend(['Model response', 'Flight measurements'])    
    
    

def compare_a(tstart, tend):

    #** set parameters **
    b0 = 0.
    psi0 = Psi[(TimeData >= tstart) & (TimeData <= tend)][0]
    r0 = r[(TimeData >= tstart) & (TimeData <= tend)][0]
    p0 = p[(TimeData >= tstart) & (TimeData <= tend)][0]
    V_0 = VTAS[(TimeData >= tstart) & (TimeData <= tend)][0]
    
    hp0 = hp[(TimeData >= tstart) & (TimeData <= tend)][0] 
    m = m0 - Fuel_used[(TimeData >= tstart) & (TimeData <= tend)][0]
    rho = rho0 * np.power( ((1+(lambdaa * hp0 / Temp0))), (-((g / (lambdaa*R)) + 1))) 
    W = m * g
    mub = m / (rho * S * b)
    CL = 2 * W / (rho * V_0 ** 2 * S)
    
    timedomain = TimeData[(TimeData >= tstart) & (TimeData <= tend)]
    timeinterval = timedomain[-1] - timedomain[0]
    
    da = d_a[(TimeData >= tstart) & (TimeData <= tend)] - d_a[(TimeData >= tstart) & (TimeData <= tend)][0]
    dr = d_r[(TimeData >= tstart) & (TimeData <= tend)] - d_r[(TimeData >= tstart) & (TimeData <= tend)][0]
    
    d_in = np.vstack((da,dr)) * np.pi/180.
    
    model = Asymmetric_Model_numerical(V_0, b, CYb, mub, CL, CYp, CYr, CYda, CYdr, Clb, KZ2, Cnb, KXZ, KX2, Clp, Cnp, Clr, Cnr, Clda, Cnda, Cldr, Cndr)
    tt_model, out_model = model.forced(d_in, timeinterval, 0.1, 0, 0, 0, 0)
    
    #**plotting**
    plt.figure()
    plt.figlegend('Model response')
    
    plt.subplot(231) #Roll angle
    plt.title('Roll angle [°]')
    plt.plot(timedomain, -out_model[1])
    plt.plot(timedomain, Psi[(TimeData >= tstart) & (TimeData <= tend)] - psi0)
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [°]')
    
    plt.subplot(232) #Roll Rate
    plt.title('Roll rate')
    plt.plot(timedomain, -out_model[2])
    plt.plot(timedomain, p[(TimeData >= tstart) & (TimeData <= tend)] - p0)
    plt.xlabel('Time [s]')
    plt.ylabel('Angular velocity [°/s]')
    
    plt.subplot(233) #Yaw rate
    plt.title('Yaw rate [°/s]')
    plt.plot(timedomain, -out_model[3])
    plt.plot(timedomain, r[(TimeData >= tstart) & (TimeData <= tend)]- r0)
    plt.xlabel('Time [s]')
    plt.ylabel('Angular velocity [°/s]')
    
    plt.subplot(234) #Aileron
    plt.title('Measured aileron deflection [°]')
    plt.plot(timedomain, da, 'darkorange')
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [°]')
    
    plt.subplot(235) #Rudder
    plt.title('Measured rudder deflection [°]')
    plt.plot(timedomain, dr, 'darkorange')
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [°]')
    
    plt.figlegend(['Model response', 'Flight measurements'])    


def compare_s_matched(tstart, tend):
    
    #** set parameters **
    a0 = VaneAOA[(TimeData >= tstart) & (TimeData <= tend)][0]
    theta0 = Pitch[(TimeData >= tstart) & (TimeData <= tend)][0]
    V_0 = VTAS[(TimeData >= tstart) & (TimeData <= tend)][0]
    q0 = q[(TimeData >= tstart) & (TimeData <= tend)][0]
    
    hp0 = hp[(TimeData >= tstart) & (TimeData <= tend)][0] 
    m = m0 - Fuel_used[(TimeData >= tstart) & (TimeData <= tend)][0]
    rho = rho0 * np.power( ((1+(Lambda * hp0 / Temp0))), (-((g / (Lambda*R)) + 1))) 
    
    W = m * g
    muc = m / (rho * S * c)
    print(V_0)
    CL = 2 * W / (rho * V_0 ** 2 * S)
    CD = CD0 + (CLa * a0) ** 2 / (np.pi * A * e)
    CX0 = CL * np.sin(theta0*np.pi/180.)
    CZ0 = -CL * np.cos(theta0*np.pi/180.)
    print(CX0, CZ0)
    q_in = q0*c/V_0
    
    timedomain = TimeData[(TimeData >= tstart) & (TimeData <= tend)]
    de = d_e[(TimeData >= tstart) & (TimeData <= tend)] - d_e[(TimeData >= tstart) & (TimeData <= tend)][0]
    det = d_et[(TimeData >= tstart) & (TimeData <= tend)]
    de_in = np.vstack((de, det)) * np.pi/180.
    
    timeinterval = timedomain[-1] - timedomain[0]

    model = Symmetrical_Model_Numerical(V_0, c, CXu, CXa, CZ0, CZu, CZa, CZadot, muc, CZq, Cmu, Cmadot, KY2, Cma, CX0, Cmq, CXde, CXdt, CZde, CZdt, Cmde, Cmdt)
    
    tt_model, out_model = model.forced(de_in, timeinterval, 0.1, 0, 0, 0, 0 , hp0)
    
    model_m = Symmetrical_Model_Numerical(V_0, c, CXu, CXa, CZ0, CZu, CZa, CZadot, muc, CZq, Cmu, Cmadot, KY2, -0.43, CX0, Cmq, CXde, CXdt, CZde, CZdt, -1.0, Cmdt)

    tt_matched, out_matched = model_m.forced(de_in, timeinterval, 0.1, 0, 0, 0, 0 , hp0)
    
    #**plotting**
    
    plt.figure()
    plt.figlegend('Model response')
    
    plt.subplot(231) #x Velocity
    plt.title('x velocity')
    plt.plot(timedomain, out_model[0], label = 'Original model response')
    plt.plot(timedomain, VTAS[(TimeData >= tstart) & (TimeData <= tend)])
    plt.plot(timedomain, out_matched[0], label = 'Corrected model response')
    plt.xlabel('Time [s]')
    plt.ylabel('Velocity [m/s]')
    
    plt.subplot(232) #AOA
    plt.title('Angle of attack')
    plt.plot(timedomain, out_model[1])
    plt.plot(timedomain, VaneAOA[(TimeData >= tstart) & (TimeData <= tend)] - a0)
    plt.plot(timedomain, out_matched[1])
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [°]')
    
    plt.subplot(233) #Pitch
    plt.title('Pitch attitude')
    plt.plot(timedomain, out_model[2])
    plt.plot(timedomain, Pitch[(TimeData >= tstart) & (TimeData <= tend)] - theta0)
    plt.plot(timedomain, out_matched[2])
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [°]')
    
    plt.subplot(234) #q
    plt.title('Pitch rate')
    plt.plot(timedomain, out_model[3])
    plt.plot(timedomain, q[(TimeData >= tstart) & (TimeData <= tend)]- q0)
    plt.plot(timedomain, out_matched[3])
    plt.xlabel('Time [s]')
    plt.ylabel('Angular velocity [°/s]')
    
    plt.subplot(235) #Altitude
    plt.title('Altitude')
    plt.plot(timedomain, out_model[4])
    plt.plot(timedomain, hp[(TimeData >= tstart) & (TimeData <= tend)])
    plt.plot(timedomain, out_matched[4])
    plt.xlabel('Time [s]')
    plt.ylabel('Altitude [m]')
    
    plt.subplot(236) #Elevator
    plt.title('Measured elevator deflection')
    plt.plot(timedomain, de, 'darkorange')
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [°]')
    
    plt.figlegend(['Original model response', 'Flight measurements', 'Corrected model response'])    
    
def compare_a_matched(tstart, tend):

    #** set parameters **
    b0 = 0.
    psi0 = Psi[(TimeData >= tstart) & (TimeData <= tend)][0]
    r0 = r[(TimeData >= tstart) & (TimeData <= tend)][0]
    p0 = p[(TimeData >= tstart) & (TimeData <= tend)][0]
    V_0 = VTAS[(TimeData >= tstart) & (TimeData <= tend)][0]
    
    hp0 = hp[(TimeData >= tstart) & (TimeData <= tend)][0] 
    m = m0 - Fuel_used[(TimeData >= tstart) & (TimeData <= tend)][0]
    rho = rho0 * np.power( ((1+(Lambda * hp0 / Temp0))), (-((g / (Lambda*R)) + 1))) 
    W = m * g
    mub = m / (rho * S * b)
    CL = 2 * W / (rho * V_0 ** 2 * S)
    
    timedomain = TimeData[(TimeData >= tstart) & (TimeData <= tend)]
    timeinterval = timedomain[-1] - timedomain[0]
    
    da = d_a[(TimeData >= tstart) & (TimeData <= tend)] - d_a[(TimeData >= tstart) & (TimeData <= tend)][0]
    dr = d_r[(TimeData >= tstart) & (TimeData <= tend)] - d_r[(TimeData >= tstart) & (TimeData <= tend)][0]
    
    d_in = np.vstack((da,dr)) * np.pi/180.
    
    model = Asymmetric_Model_numerical(V_0, b, CYb, mub, CL, CYp, CYr, CYda, CYdr, Clb, KZ2, Cnb, KXZ, KX2, Clp, Cnp, Clr, Cnr, Clda, Cnda, Cldr, Cndr)
    tt_model, out_model = model.forced(d_in, timeinterval, 0.1, 0, 0, 0, 0)
    
    model_matched = Asymmetric_Model_numerical(V_0, b, CYb, mub, CL, CYp, CYr, CYda, CYdr, -0.12260, KZ2, +0.1048, KXZ, KX2, Clp, Cnp, Clr, Cnr, Clda, -0.0470, Cldr, -0.0739)
    tt_matched, out_matched = model_matched.forced(d_in, timeinterval, 0.1, 0, 0, 0, 0)
    
    #**plotting**
    plt.figure()
    
    plt.subplot(231) #Roll angle
    plt.title('Roll angle [°]')
    plt.plot(timedomain, -out_model[1])
    plt.plot(timedomain, Psi[(TimeData >= tstart) & (TimeData <= tend)] - psi0)
    plt.plot(timedomain, -out_matched[1])
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [°]')
    
    plt.subplot(232) #Roll Rate
    plt.title('Roll rate')
    plt.plot(timedomain, -out_model[2])
    plt.plot(timedomain, p[(TimeData >= tstart) & (TimeData <= tend)] - p0)
    plt.plot(timedomain, -out_matched[2])
    plt.xlabel('Time [s]')
    plt.ylabel('Angular velocity [°/s]')
    
    plt.subplot(233) #Yaw rate
    plt.title('Yaw rate [°/s]')
    plt.plot(timedomain, -out_model[3])
    plt.plot(timedomain, r[(TimeData >= tstart) & (TimeData <= tend)]- r0)
    plt.plot(timedomain, -out_matched[3])
    plt.xlabel('Time [s]')
    plt.ylabel('Angular velocity [°/s]')
    
    plt.subplot(234) #Aileron
    plt.title('Measured aileron deflection [°]')
    plt.plot(timedomain, da, 'darkorange')
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [°]')
    
    plt.subplot(235) #Rudder
    plt.title('Measured rudder deflection [°]')
    plt.plot(timedomain, dr, 'darkorange')
    plt.xlabel('Time [s]')
    plt.ylabel('Angle [°]')
    
    plt.figlegend(['Original model response', 'Flight measurements', 'Corrected model response'])       
    

    
    
compare_s_matched(2805, 2950)
compare_s_matched(3100, 3125)
compare_s(3020,3050)
compare_a_matched(3180, 3250)
compare_a_matched(3425, 3600)
compare_a_matched(3020, 3050)

   