# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:14:25 2019

@author: M.Osowski
"""

import numpy as np
import matplotlib.pyplot as plt
import control 

from Cit_par import *

V = V0

class Symmetrical_Model_Numerical():
    def __init__(self):
        
        #Setup of A matrix
        xu = V * CXu / (c * 2 * muc);
        xa = V * CXa / (c * 2 * muc);
        x0 = V * CZ0 / (c * 2 * muc);
            
        zu = V * CZu / (c * (2*muc - CZadot));
        za = V * CZa / (c * (2*muc - CZadot));
        z0 = V * CZ0 / (c * (2*muc - CZadot));
        zq = V * CZq / (c * (2*muc - CZadot));
            
        mu = V * (Cmu + CZu * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2);
        ma = V * (Cma + CZa * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2);
        m0 = V * (CX0 * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2);
        mq = V * (Cmq + Cmadot * (2 * muc + CZq) / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
         
        self.A = np.array(([0, V, 0, 0 ,0], [0, xu, xa, x0, 0],[0, zu, za, z0, zq], [0, 0, 0, 0, V0/c], [0, mu, ma, m0, mq]))
    
        #Setup of B matrix
        xde = V * CXde / (c * 2 * muc)
        xdt = V * CXdt / (c * 2 * muc)
        
        zde = V * CZde / (c * (2 * muc - CZadot))
        zdt = V * CZdt / (c * (2 * muc - CZadot))
        
        mde = V * (Cmde + CZde * (2 * Cmadot) / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
        mdt = V * (Cmdt + CZdt * (2 * Cmadot) / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
        
        self.B = np.array(([0, 0], [xde, xdt], [zde, zdt], [0, 0], [mde, mdt]))
        
        self.C = np.array(([1, 0, 0, 0, 0], [0, V, 0, 0, 0], [0, 0, 1, 0, 0]))
        
        self.D = np.array(([0, 0], [0, 0], [0, 0]))
        
        self.sys = control.ss(self.A, self.B, self.C, self.D)
        
    def pulse_e(self, d_e = 10. * np.pi / 180., tmax = 30., step = 0.1):
        tt = np.arange(0, tmax, step)
        
        control.impulse_response(self.sys, tt, d_e, 0, (0,1,2))


citation = Symmetrical_Model_Numerical()

print(citation.sys)
