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


def plot_data(tstart, tend, datasets, a0, theta0):
    
    Flightdata = np.io.loadmat('FTISxprt-20190319_123022.mat')

    TimeData = Flightdata['flightdata']['time'][0][0][0][0][0].transpose()
    VaneAOA = Flightdata['flightdata']['vane_AOA'][0][0][0][0][0] - a0
    Pitch = Flightdata['flightdata']['Ahrs1_Pitch'][0][0][0][0][0] - theta0
    VTAS = Flightdata['flightdata']['Dadc1_tas'][0][0][0][0][0] * 0.514444 #converted to m/s
       
    
    timedomain = TimeData[(TimeData >= tstart) & (TimeData <= tend)] - tstart
    
    if type(datasets) == list:
        for Set in datasets:
            measurements = Set[(TimeData >= tstart) & (TimeData <= tend)]
            plt.plot(timedomain, measurements)
    else:
        measurements = datasets[(TimeData >= tstart) & (TimeData <= tend)]
        plt.plot(timedomain, measurements)
        
    plt.show()
    

#plot_data(2790,3000,Pitch, 0 ,0)




