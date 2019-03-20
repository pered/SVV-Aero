#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:30:16 2019

@author: tristan
"""

import control
import matplotlib.pyplot as plt
from Cit_par import *

V = V0

class Asymmetric_Model_numerical():
    
    def __init__(self):
        
        self.set_A()
        self.set_B()
        self.set_C()
        self.set_D()
    
    y_beta   = (V/b) * (CYb/(2*mub))
    y_phi    = (V/b) * (CL/(2*mub))
    y_p      = (V/b) * (CYp/(2*mub))
    y_r      = (V/b) * ((CYr-4*mub)/(2*mub))
    y_deltaa = (V/b) * (CYda/(2*mub))
    y_deltar = (V/b) * (CYdr/(2*mub))
    
    l_beta   = (V/b) * ((Clb * KZ2 + Cnb * KXZ)/(4*mub * (KX2 * KZ2 - KXZ**2)))
    l_phi    = 0
    l_p      = (V/b) * ((Clp * KZ2 + Cnp * KXZ)/(4*mub * (KX2 * KZ2 - KXZ**2)))
    l_r      = (V/b) * ((Clr * KZ2 + Cnr * KXZ)/(4*mub * (KX2 * KZ2 - KXZ**2)))
    l_deltaa = (V/b) * ((Clda * KZ2 + Cnda * KXZ)/(4*mub * (KX2 * KZ2 - KXZ**2)))
    l_deltar = (V/b) * ((Cldr * KZ2 + Cndr * KXZ)/(4*mub * (KX2 * KZ2 - KXZ**2)))
    
    n_beta   = (V/b) * ((Clb * KXZ + Cnb * KX2)/(4*mub * (KX2 * KZ2 - KXZ**2)))
    n_phi    = 0
    n_p      = (V/b) * ((Clp * KXZ + Cnp * KX2)/(4*mub * (KX2 * KZ2 - KXZ**2)))
    n_r      = (V/b) * ((Clr * KXZ + Cnr * KX2)/(4*mub * (KX2 * KZ2 - KXZ**2)))
    n_deltaa = (V/b) * ((Clda * KXZ + Cnda * KX2)/(4*mub * (KX2 * KZ2 - KXZ**2)))
    n_deltar = (V/b) * ((Cldr * KXZ + Cndr * KX2)/(4*mub * (KX2 * KZ2 - KXZ**2)))

    def set_A(self):
        
        self.A = [[ self.y_beta, self.y_phi,       self.y_p, self.y_r ],
                  [          0.,         0.,      2 * (V/b),       0. ],
                  [ self.l_beta,         0.,       self.l_p, self.l_r ],
                  [ self.n_beta,         0.,       self.n_p, self.n_r ]]
        
    def set_B(self):
        
        self.B = [[            0., self.y_deltar],
                  [            0.,            0.],
                  [ self.l_deltaa, self.l_deltar],
                  [ self.n_deltaa, self.n_deltar]]
        
    def set_C(self):
        
        self.C = [[  1., 0.,      0.,      0. ],
                  [  0., 1.,      0.,      0. ],
                  [  0., 0., (2*V)/b,      0. ],
                  [  0., 0.,      0., (2*V)/b ],
                  [ -1., 0.,      0.,      0. ]]
        
    def set_D(self):
        
        self.D = [[ 0., 0. ],
                  [ 0., 0. ],
                  [ 0., 0. ],
                  [ 0., 0. ],
                  [ 0., 0. ]]
        
asymm = Asymmetric_Model_numerical()
        
sys1 = control.ss(asymm.A,asymm.B,asymm.C,asymm.D)

inpt = [[ 0.025],
        [ 0.]]

response = control.impulse_response(sys1,input=1)

plt.plot(response[1][5])