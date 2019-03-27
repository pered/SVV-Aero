# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 09:14:25 2019

@author: M.Osowski
"""

import scipy as np
import matplotlib.pyplot as plt
import control 
#from Cit_par import *


#plt.close('all')

#The Following class defines a model of the longitudual dynamics of an aircraft. A state space model is used, with
#state variables u^, alpha, theta, and qc/V. The outputs of the state space system are u, alpha and theta. Parameters
#from the imported parameter file are used to construct the state space and provide initial conditions. The model 
#contains several input modes, which can be called as class functions. These plot the response of the aircraft to a
#given disturbance

class Symmetrical_Model_Numerical():
    def __init__(self, V, c, CXu, CXa, CZ0, CZu, CZa, CZadot, muc, CZq, Cmu, Cmadot, KY2, Cma, CX0, Cmq, CXde, CXdt, CZde, CZdt, Cmde, Cmdt):
        #Upon creation of the class object, the state space system is set up and assembled
        self.V = V
        self.c = c
        #Setup of A matrix
        self.xu = V * CXu / (c * 2 * muc)
        self.xa = V * CXa / (c * 2 * muc)
        self.x0 = V * CZ0 / (c * 2 * muc)
            
        self.zu = V * CZu / (c * (2*muc - CZadot))
        self.za = V * CZa / (c * (2*muc - CZadot))
        self.z0 = -V * CX0 / (c * (2*muc - CZadot))
        self.zq = V * (2*muc + CZq) / (c * (2*muc - CZadot))
            
        self.mu = V * (Cmu + CZu * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
        self.ma = V * (Cma + CZa * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
        self.m0 = -V * (CX0 * Cmadot / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
        self.mq = V * (Cmq + Cmadot * (2 * muc + CZq) / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
         
        self.A = np.array(([self.xu, self.xa, self.x0, 0, 0],
                           [self.zu, self.za, self.z0, self.zq, 0],
                           [0, 0, 0, V/c, 0], 
                           [self.mu, self.ma, self.m0, self.mq, 0],
                           [0, 0, V, 0, 0]))
    
        #Setup of B matrix
        xde = V * CXde / (c * 2 * muc)
        xdt = V * CXdt / (c * 2 * muc)
        
        zde = V * CZde / (c * (2 * muc - CZadot))
        zdt = V * CZdt / (c * (2 * muc - CZadot))
        
        mde = V * (Cmde + CZde * (Cmadot) / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
        mdt = V * (Cmdt + CZdt * (Cmadot) / (2 * muc - CZadot)) / (c * 2 * muc * KY2)
        
        self.B = np.array(([xde, xdt], 
                           [zde, zdt], 
                           [0, 0], 
                           [mde, mdt],
                           [0, 0]))
        
        #Setup of C matrix
        print(V,c)
        self.C = np.array(([V, 0, 0, 0, 0],
                           [0, 1, 0, 0, 0],
                           [0, 0, 1, 0, 0],
                           [0, 0, 0, V/c, 0],
                           [0, 0, 0, 0, 1]))
        
        #Setup of D matrix
        self.D = np.array(([0, 0], 
                           [0, 0], 
                           [0, 0],
                           [0, 0],
                           [0, 0]))
        
        #Assembly of state space system
        self.sys = control.ss(self.A, self.B, self.C, self.D)
    
    def ic(self, tmax, step, u0, a0, theta0, q0, hp0):
        tt = np.arange(0, tmax + step, step)
        x0 = np.array([u0, a0 * np.pi / 180., theta0 * np.pi / 180., q0 * self.V / self.c, hp0])
        return tt, x0
        
    def pulse_e(self, d_e = 10. * np.pi / 180., tmax = 180., step = 0.01, u0 = 0., a0 = 0., theta0 = 0., q0 = 0, mag = 0.04):
        
        tt, x0 = Symmetrical_Model_Numerical.ic(self, tmax, step, u0, a0, theta0, q0, hp0)
        
        tt, out = control.impulse_response(self.sys, tt, x0, 0)
        
        out *= mag
        
        V0vec = self.V * np.ones(len(tt))
        out[0] += V0vec
        out[1] *= 180/np.pi 
        out[2] *= 180/np.pi
        out[3] *= 180/np.pi
        
        print(out[1])
        
        Symmetrical_Model_Numerical.plotting(self, tt, out)
        
    def initial(self, tmax = 300., step = 0.1, u0 = 0.2, a0 = 0., theta0 = 0., q0 = 0):
        
        tt, x0 = Symmetrical_Model_Numerical.ic(self, tmax, step, u0, a0, theta0, q0, hp0)
        
        tt, out = control.initial_response(self.sys, tt, x0)
        
        V0vec = self.V * np.ones(len(tt))
        out[0] += V0vec
        out[1] *= 180/np.pi
        out[2] *= 180/np.pi
        out[3] *= 180/np.pi
        
        Symmetrical_Model_Numerical.plotting(self, tt, out)
        
        return tt, out
    
    def step(self, tmax = 300., step = 0.1, u0 = 0., a0 = 0., theta0 = 0., q0 = 0, mag = 0.02):
        
        tt, x0 = Symmetrical_Model_Numerical.ic(self, tmax, step, u0, a0, theta0, q0, hp0)
        
        tt, out = control.step_response(self.sys, tt, x0)
        
        out *= mag
        
        V0vec = self.V * np.ones(len(tt))
        out[0] += V0vec
        out[1] *= 180/np.pi
        out[2] *= 180/np.pi
        out[3] *= 180/np.pi
        
        Symmetrical_Model_Numerical.plotting(self, tt, out)
        
        return tt, out

    def pulse(self, tmax = 1000., step = 0.1, u0 = 0., a0 = 0., theta0 = 0., q0 = 0, mag = -0.02, length = 1000.):
        
        tt, x0 = Symmetrical_Model_Numerical.ic(self, tmax, step, u0, a0, theta0, q0, hp0)
        
        f = np.zeros((2, len(tt)))
        f[0, tt<=length] = mag
        
        tt, out, xout = control.forced_response(self.sys, tt, f, x0)
        
        V0vec = self.V * np.ones(len(tt))
        out[0] += V0vec
        out[1] *= 180/np.pi
        out[2] *= 180/np.pi
        out[3] *= 180/np.pi
        
        Symmetrical_Model_Numerical.plotting(self, tt, out)
        
        return tt, out
        
    def forced(self, f, tmax = 1800., step = 0.1, u0 = 0., a0 = 10., theta0 = 0., q0 = 0, hp0 = 0, plot = False):
        tt, x0 = Symmetrical_Model_Numerical.ic(self, tmax, step, u0, a0, theta0, q0, hp0)
        
        
        tt, out, xout = control.forced_response(self.sys, tt, f, x0)
        V0vec = self.V * np.ones(len(tt))
        out[0] += V0vec
        out[1] *= 180/np.pi
        out[2] *= 180/np.pi
        out[3] *= 180/np.pi
        
        
        
        if plot:
            Symmetrical_Model_Numerical.plotting(self, tt, out)
    #       
        
        return tt, out
    
    def plotting(self, tt, output):
        
        plt.figure()
        
        plt.subplot(231) #x Velocity
        plt.title('x velocity')
        plt.plot(tt, output[0], label = 'Model response')
        plt.xlabel('Time [s]')
        plt.ylabel('Velocity [m/s]')
        
        plt.subplot(232) #AOA
        plt.title('Angle of attack')
        plt.plot(tt, output[1])
        plt.xlabel('Time [s]')
        plt.ylabel('Angle [°]')
        
        plt.subplot(233) #Pitch
        plt.title('Pitch attitude')
        plt.plot(tt, output[2])
        plt.xlabel('Time [s]')
        plt.ylabel('Angle [°]')
        
        plt.subplot(234) #q
        plt.title('Pitch rate')
        plt.plot(tt, output[3])
        plt.xlabel('Time [s]')
        plt.ylabel('Angular velocity [°/s]')
        
        plt.subplot(235) #Altitude
        plt.title('Altitude')
        plt.plot(tt, output[4])
        plt.xlabel('Time [s]')
        plt.ylabel('Altitude [m]')
        
        
        plt.figlegend(['Model response'])    
        
    def eigs(self):
        self.eigs = np.linalg.eigvals(self.A)
        
        X = [x.real for x in self.eigs]
        Y = [x.imag for x in self.eigs]
        
        plt.figure()
        plt.scatter(X,Y)
        plt.xlabel('Re')
        plt.ylabel('Im')
        plt.grid()
        plt.show()
        
        print(self.eigs)
        
class Asymmetric_Model_numerical():
    
    def __init__(self, V, b, CYb, mub, CL, CYp, CYr, CYda, CYdr, Clb, KZ2, Cnb, KXZ, KX2, Clp, Cnp, Clr, Cnr, Clda, Cnda, Cldr, Cndr):
        
        self.V = V
        self.b = b
        
        self.y_beta   = (V/b) * (CYb/(2*mub))
        self.y_phi    = (V/b) * (CL/(2*mub))
        self.y_p      = (V/b) * (CYp/(2*mub))
        self.y_r      = (V/b) * ((CYr-4*mub)/(2*mub))
        self.y_deltaa = (V/b) * (CYda/(2*mub))
        self.y_deltar = (V/b) * (CYdr/(2*mub))
        
        
        self.l_beta   = (V/b) * ((Clb * KZ2 + Cnb * KXZ)/(4*mub * (KX2 * KZ2 - KXZ**2)))
        self.l_phi    = 0
        self.l_p      = (V/b) * ((Clp * KZ2 + Cnp * KXZ)/(4*mub * (KX2 * KZ2 - KXZ**2)))
        self.l_r      = (V/b) * ((Clr * KZ2 + Cnr * KXZ)/(4*mub * (KX2 * KZ2 - KXZ**2)))
        self.l_deltaa = (V/b) * ((Clda * KZ2 + Cnda * KXZ)/(4*mub * (KX2 * KZ2 - KXZ**2)))
        self.l_deltar = (V/b) * ((Cldr * KZ2 + Cndr * KXZ)/(4*mub * (KX2 * KZ2 - KXZ**2)))
        
        self.n_beta   = (V/b) * ((Clb * KXZ + Cnb * KX2)/(4*mub * (KX2 * KZ2 - KXZ**2)))
        self.n_phi    = 0
        self.n_p      = (V/b) * ((Clp * KXZ + Cnp * KX2)/(4*mub * (KX2 * KZ2 - KXZ**2)))
        self.n_r      = (V/b) * ((Clr * KXZ + Cnr * KX2)/(4*mub * (KX2 * KZ2 - KXZ**2)))
        self.n_deltaa = (V/b) * ((Clda * KXZ + Cnda * KX2)/(4*mub * (KX2 * KZ2 - KXZ**2)))
        self.n_deltar = (V/b) * ((Cldr * KXZ + Cndr * KX2)/(4*mub * (KX2 * KZ2 - KXZ**2)))
        
        self.set_A()
        self.set_B()
        self.set_C()
        self.set_D()
        
        self.sys = control.ss(self.A,self.B,self.C,self.D)

    def set_A(self):
        
        self.A = [[ self.y_beta, self.y_phi,       self.y_p, self.y_r ],
                  [          0.,         0.,      2 * (self.V/self.b),       0. ],
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
                  [  0., 0., (2*self.V)/self.b,      0. ],
                  [  0., 0.,      0., (2*self.V)/self.b ],
                  [ -1., 0.,      0.,      0. ]]
        
    def set_D(self):
        
        self.D = [[ 0., 0. ],
                  [ 0., 0. ],
                  [ 0., 0. ],
                  [ 0., 0. ],
                  [ 0., 0. ]] 
        
    def ic(self, tmax, step, b0, psi0, p0, r0):
        tt = np.arange(0, tmax + step, step)
        x0 = np.array([b0 * np.pi / 180., psi0 * np.pi / 180., p0*self.b/(2*self.V), r0 *self.b/(2*self.V)])
        
        return tt, x0
    
    def initial(self, tmax = 300., step = 0.1, b0 = 5., psi0 = 0., r0 = 0, p0 = 0):
        
        tt, x0 = Asymmetric_Model_numerical.ic(self, tmax, step, b0, psi0, p0, r0)
        tt, out = control.initial_response(self.sys, tt, x0)
        out *= 180/np.pi
        Asymmetric_Model_numerical.plotting(tt, out)
    
    def forced(self, f, tmax = 300., step = 0.1, b0 = 0, psi0 = 0, r0 = 0, p0 = 0, plot = False):
        
        tt, x0 = Asymmetric_Model_numerical.ic(self, tmax, step, b0, psi0, p0, r0)
        tt, out, xout = control.forced_response(self.sys, tt, f, x0)
        out *= 180/np.pi
        if plot:    
            Asymmetric_Model_numerical.plotting(tt, out)
        
        return tt, out
    
    def pulse(self, tmax = 100., step = 0.1, b0 = 0, psi0 = 0, r0 = 0, p0 = 0, mag = -0.02, length = 100., defl = 'aileron'):
        
        tt, x0 = Asymmetric_Model_numerical.ic(self, tmax, step, b0, psi0, p0, r0)
        f = np.zeros((2, len(tt)))
        if defl == 'rudder':
            f[1, tt<=length] = mag
        if defl == 'aileron':
            f[0, tt<=length] = mag
        tt, out, xout = control.forced_response(self.sys, tt, f, x0)
        out *= 180/np.pi
        Asymmetric_Model_numerical.plotting(tt, out)
        
        return tt, out
    
    def plotting(tt, output):
        
        plt.figure()
        
        plt.subplot(221) #Roll angle
        plt.title('Roll angle')
        plt.plot(tt, -output[1])
        plt.xlabel('Time [s]')
        plt.ylabel('Angle [°]')
        
        plt.subplot(222) #Roll Rate
        plt.title('Roll rate')
        plt.plot(tt, -output[2])
        plt.xlabel('Time [s]')
        plt.ylabel('Angular velocity [°/s]')
        
        plt.subplot(223) #Sideslip
        plt.title('Sideslip')
        plt.plot(tt, -output[0])
        plt.xlabel('Time [s]')
        plt.ylabel('Angle [°]')
        
        plt.subplot(224) #Yaw rate
        plt.title('Yaw rate')
        plt.plot(tt, -output[3])
        plt.xlabel('Time [s]')
        plt.ylabel('Angular velocity [°/s]')
        
        plt.figlegend(['Model response'])    

        

#citation = Symmetrical_Model_Numerical(V0, c, CXu, CXa, CZ0, CZu, CZa, CZadot, muc, CZq, Cmu, Cmadot, KY2, Cma, CX0, Cmq, CXde, CXdt, CZde, CZdt, Cmde, Cmdt)
#citation.pulse_e()
#citation.initial()
#citation.step()
#citation.pulse()
#citation.eigs()
        
#citation_a = Asymmetric_Model_numerical(V0, b, CYb, mub, CL, CYp, CYr, CYda, CYdr, Clb, KZ2, Cnb, KXZ, KX2, Clp, Cnp, Clr, Cnr, Clda, Cnda, Cldr, Cndr)
#citation_a.pulse()
#citation_a.initial()