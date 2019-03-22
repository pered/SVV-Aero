# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 15:27:49 2019

@author: M.Osowski
"""

from Cit_par_refdata import *
import scipy as np
import matplotlib.pyplot as plt
import control
from Cit_class_experimental import Symmetrical_Model_Numerical
from Flight_data_reader import plot_data


 
    
Flightdata = np.io.loadmat('FTISxprt-20190319_123022.mat')

TimeData = Flightdata['flightdata']['time'][0][0][0][0][0].transpose()
VaneAOA = Flightdata['flightdata']['vane_AOA'][0][0][0][0][0]
Pitch = Flightdata['flightdata']['Ahrs1_Pitch'][0][0][0][0][0]
d_e = Flightdata['flightdata']['delta_e'][0][0][0][0][0]

def compare(tstart, tend):
    timedomain = TimeData[(TimeData >= tstart) & (TimeData <= tend)]
    de = d_e[(TimeData >= tstart) & (TimeData <= tend)] 
    de_in = np.vstack((de, np.zeros(len(de)))) * np.pi/180
    a0 = VaneAOA[(TimeData >= tstart) & (TimeData <= tend)][0]
    theta0 = Pitch[(TimeData >= tstart) & (TimeData <= tend)][0]
    timeinterval = timedomain[-1] - timedomain[0]
    
    model = Symmetrical_Model_Numerical(V0)
    
    model.forced(de_in, timeinterval, 0.1, 0, a0, theta0, 0)
    
    plt.plot(timedomain - tstart, de)
    
    plot_data(tstart, tend, [Pitch])


    




compare(200, 400)


