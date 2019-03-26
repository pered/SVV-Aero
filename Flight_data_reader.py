# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 14:03:35 2019

@author: M.Osowski
"""

import scipy as np
import matplotlib.pyplot as plt

#Flightdata = np.io.loadmat('FTISxprt-20190319_123022.mat')
#
#TimeData = Flightdata['flightdata']['time'][0][0][0][0][0].transpose()
#VaneAOA = Flightdata['flightdata']['vane_AOA'][0][0][0][0][0]
#Pitch = Flightdata['flightdata']['Ahrs1_Pitch'][0][0][0][0][0]


def plot_data(tstart, tend, datasets):
    
    Flightdata = np.io.loadmat('FTISxprt-20190319_123022.mat')

    TimeData = Flightdata['flightdata']['time'][0][0][0][0][0].transpose()
    VaneAOA = Flightdata['flightdata']['vane_AOA'][0][0][0][0][0]
    Pitch = Flightdata['flightdata']['Ahrs1_Pitch'][0][0][0][0][0]
    VTAS = Flightdata['flightdata']['Dadc1_tas'][0][0][0][0][0] * 0.514444 #converted to m/s
    Psi = Flightdata['flightdata']['Ahrs1_Roll'][0][0][0][0][0]
    p = Flightdata['flightdata']['Ahrs1_bRollRate'][0][0][0][0][0]
    r = Flightdata['flightdata']['Ahrs1_bYawRate'][0][0][0][0][0]
    
    d_e = Flightdata['flightdata']['delta_e'][0][0][0][0][0]
    d_et = Flightdata['flightdata']['elevator_dte'][0][0][0][0][0]
    d_ec = d_e + d_et
       
    
    timedomain = TimeData[(TimeData >= tstart) & (TimeData <= tend)] - tstart
    
    if type(datasets) == list:
        for Set in datasets:
            measurements = Set[(TimeData >= tstart) & (TimeData <= tend)]
            plt.plot(timedomain, measurements)
    else:
        measurements = datasets[(TimeData >= tstart) & (TimeData <= tend)]
        plt.plot(timedomain, measurements)
        
    plt.show()
    
    

#plot_data(0,45000,[Pitch,d_e])




