# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:14:25 2019

@author: M.Osowski
"""

import scipy as np
import matplotlib.pyplot as plt
import control 

from Cit_par_refdata import *



plt.close('all')

#The Following class defines a model of the longitudual dynamics of an aircraft. A state space model is used, with
#state variables u^, alpha, theta, and qc/V. The outputs of the state space system are u, alpha and theta. Parameters
#from the imported parameter file are used to construct the state space and provide initial conditions. The model 
#contains several input modes, which can be called as class functions. These plot the response of the aircraft to a
#given disturbance

class Symmetrical_Model_Numerical():
    def __init__(self, V):
        #Upon creation of the class object, the state space system is set up and assembled
        self.V = V
        #Setup of A matrix
        self.xu = V * CXu / (c * 2 * muc)
        self.xa = V * CXa / (c * 2 * muc)
        self.x0 = V * CZ0 / (c * 2 * muc)
            
        self.zu = V * CZu / (c * (2*muc - CZadot))
        self.za = V * CZa / (c * (2*muc - CZadot))
        self.z0 = -V * CZ0 / (c * (2*muc - CZadot))
        self.zq = V * (2*muc + CZq) / (c * (2*muc - CZadot))
            
        self.mu = V * (Cmu + CZu * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
        self.ma = V * (Cma + CZa * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
        self.m0 = -V * (CX0 * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
        self.mq = V * (Cmq + Cmadot * (2 * muc + CZq) / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
         
        self.A = np.array(([self.xu, self.xa, self.x0, 0],[self.zu, self.za, self.z0, self.zq], [0, 0, 0, V/c], [self.mu, self.ma, self.m0, self.mq]))
    
        #Setup of B matrix
        xde = V * CXde / (c * 2 * muc)
        xdt = V * CXdt / (c * 2 * muc)
        
        zde = V * CZde / (c * (2 * muc - CZadot))
        zdt = V * CZdt / (c * (2 * muc - CZadot))
        
        mde = V * (Cmde + CZde * (2 * Cmadot) / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
        mdt = V * (Cmdt + CZdt * (2 * Cmadot) / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
        
        self.B = np.array(([xde, xdt], [zde, zdt], [0, 0], [mde, mdt]))
        
        #Setup of C matrix
        self.C = np.array(([V, 0, 0, 0], [0, 1, 0, 0],[0, 0, 1, 0]))
        
        #Setup of D matrix
        self.D = np.array(([0, 0], [0, 0], [0, 0]))
        
        #Assembly of state space system
        self.sys = control.ss(self.A, self.B, self.C, self.D)
    
    def ic(self, tmax, step, u0, a0, theta0, q0):
        tt = np.arange(0, tmax + step, step)
        x0 = np.array([u0, a0 * np.pi / 180., theta0 * np.pi / 180., q0 * self.V / c])
        
        return tt, x0
        
    def pulse_e(self, d_e = 10. * np.pi / 180., tmax = 180., step = 0.01, u0 = 0., a0 = 2., theta0 = 2., q0 = 0, mag = 0.005):
        
        tt, x0 = Symmetrical_Model_Numerical.ic(self, tmax, step, u0, a0, theta0, q0)
        
        tt, out = control.impulse_response(self.sys, tt, x0, 0)
        
        out *= mag
        
        V0vec = self.V * np.ones(len(tt))
        out[0] += V0vec
        out[1] *= 180/np.pi
        out[2] *= 180/np.pi
    
        
        #plt.close('all')   
        plt.plot(tt,out[0])
        plt.plot(tt,out[1])
        plt.plot(tt,out[2])
        plt.show()
        
    def initial(self, tmax = 180., step = 0.1, u0 = 0., a0 = 10., theta0 = 10., q0 = 0):
        
        tt, x0 = Symmetrical_Model_Numerical.ic(self, tmax, step, u0, a0, theta0, q0)
        
        tt, out = control.initial_response(self.sys, tt, x0)
        
        V0vec = self.V * np.ones(len(tt))
        out[0] += V0vec
        out[1] *= 180/np.pi
        out[2] *= 180/np.pi
        
        #plt.close('all')   
        plt.plot(tt,out[0])
        plt.plot(tt,out[1])
        plt.plot(tt,out[2])
        plt.show()
        
        return tt, out
    
    def step(self, tmax = 1800., step = 0.1, u0 = 0., a0 = 10., theta0 = 10., q0 = 0, mag = 0.01):
        
        tt, x0 = Symmetrical_Model_Numerical.ic(self, tmax, step, u0, a0, theta0, q0)
        
        tt, out = control.step_response(self.sys, tt, x0)
        
        out *= mag
        
        V0vec = self.V * np.ones(len(tt))
        out[0] += V0vec
        out[1] *= 180/np.pi
        out[2] *= 180/np.pi
        
        #plt.close('all')   
        plt.plot(tt,out[0])
        plt.plot(tt,out[1])
        plt.plot(tt,out[2])
        plt.show()
        
        return tt, out

    def pulse(self, tmax = 1800., step = 0.1, u0 = 0., a0 = 10., theta0 = 10., q0 = 0, mag = 0.01, length = 2.):
        
        tt, x0 = Symmetrical_Model_Numerical.ic(self, tmax, step, u0, a0, theta0, q0)
        
        f = np.zeros((2, len(tt)))
        f[0, tt<=length] = mag
        
        tt, out, xout = control.forced_response(self.sys, tt, f, x0)
        
        V0vec = self.V * np.ones(len(tt))
        out[0] += V0vec
        out[1] *= 180/np.pi
        out[2] *= 180/np.pi
        
        #plt.close('all')   
        plt.plot(tt,out[0])
        plt.plot(tt,out[1])
        plt.plot(tt,out[2])
        plt.show()
        
        return tt, out
        
    def forced(self, f, tmax = 1800., step = 0.1, u0 = 0., a0 = 10., theta0 = 10., q0 = 0):
        tt, x0 = Symmetrical_Model_Numerical.ic(self, tmax, step, u0, a0, theta0, q0)
        
        
        tt, out, xout = control.forced_response(self.sys, tt, f, x0)
        
        V0vec = V0 * np.ones(len(tt))
        out[0] += V0vec
        out[1] *= 180/np.pi
        out[2] *= 180/np.pi
        
        #plt.close('all')   
        plt.plot(tt,out[0])
        #plt.plot(tt,out[1])
        plt.plot(tt,out[2])
        plt.show()
        
        return tt, out
    
    def eigs(self):
        self.eigs = np.linalg.eigvals(self.A)
        
        X = [x.real for x in self.eigs]
        Y = [x.imag for x in self.eigs]
        
        plt.figure()
        plt.scatter(X,Y)
        plt.show()
        
        print(self.eigs)
        
    
#citation = Symmetrical_Model_Numerical(V0)
#citation.pulse_e()
#citation.initial()
#citation.step()
#citation.pulse()
#citation.eigs()
