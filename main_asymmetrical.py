#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:30:16 2019

@author: tristan
"""

import control
import matplotlib.pyplot as plt
from Cit_par import *
from Flight_data_reader import *

V = 59.9

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
        
"""---------------------------------------------------------"""
        
def inputs():
    
    print("Define Initial Conditions:")
    b0   = float(input("Angle of sideslip $\\beta$: "))
    phi0 = float(input("Angle of roll \\phi$: "))
    V    = float(input("Magnitude of the airspeed vector V: "))
    p0   = float(input("Angular velocity about X-axis p: "))
    r0   = float(input("Angular velocity about Z-axis r: "))

    x0 = [[             b0 ],
          [           phi0 ],
          [   (p0*b)/(2.*V) ],
          [   (r0*b)/(2.*V) ]]
    
    print("Define timeframe")
    t_end = int(input("Time t: "))
    dt = float(input("Stepsize in seconds dt: "))
    
    t = np.arange(0, int(t_end), dt)
       
    print("Define Inputs:")
    delta_a   = float(input("Aileron deflection $\\delta$_a: "))
    delta_r = float(input("Rudder deflection \\delta$_r: "))
    
    inpt = [ np.concatenate((np.full((int(1/dt)), delta_a), np.zeros(int((t_end-1)*(1/dt)))), axis=None ) ,
             np.concatenate((np.full((int(1/dt)), delta_r), np.zeros(int((t_end-1)*(1/dt)))), axis=None ) ]
    
    return x0, inpt, t

"""---------------------------------------------------------"""

def gen_plots(resp,title):   
    
    labels = ["Angle of sideslip $\\beta$ [Rad]", "Angle of roll $\\phi$ [Rad]", "Angular velocity about the X-axis p [Rad/sec]", "Angular velocity about the Z-axis r [Rad/sec]", "Yaw angle $\\psi$ [Rad]"]
    
    plt.suptitle(title)
    
    for i in range(5):
        
        plt.subplot(2,3,i+1)
        plt.plot(t, resp[1][i])
        plt.xlabel("time[s]")
        plt.ylabel(labels[i])
        
    plt.show()

"""---------------------------------------------------------"""

x0, inpt, t = inputs()

asymm = Asymmetric_Model_numerical()
        
sys1 = control.ss(asymm.A,asymm.B,asymm.C,asymm.D)

p = control.pole(sys1)
print(p)
ee = np.linalg.eig(asymm.A)

response_a = control.forced_response(sys1, U=inpt, X0=x0, T=t)
response_r = control.forced_response(sys1, U=inpt, X0=x0, T=t)


#gen_plots(response_a, "Response curves for rudder deflection")
#gen_plots(response_r, "Response curves for aileron deflection")

plot_data(46*60,48*60,response_a)