#!/usr/bin/env python

import numpy as np
from scipy.integrate import solve_ivp as ivp
import matplotlib.pyplot as plt
import Constants


# Simple Harmonic Pendulum
class SimpleHarmonicPendulum:
    # Constructor
    # length: length of pendulum [rad]
    def __init__(self, length = 1.0):
        self.len = length


    # Get simple harmonic oscillator acceleration
    # x: displacement [rad]
    # return: acceleration [rad/s]
    def get_acceleration(self, x):
        return Constants.g_acceleration * self.len * x


    # Equation of state for the SHO
    # t: time [seconds]
    # x: [displacement [rad], velocity [rad/s]]
    # return: [velocity [rad/s], acceleration [rad/s/s]]
    def deriv(self, t, x):
        return [x[1], self.get_acceleration(x[0])]


    # Solve equations of state numerically
    # y0: initial state [displacement [rad], velocity [rad/s]]
    # t:  time [s]
    # return: [time [s], displacement [rad], velocity [rad/s], acceleration [rad/s/s]]
    def solve(self, y0, t):
        sol = ivp(self.deriv, [t[0], t[-1]], y0, t_eval=t)
        return [sol.t, sol.y[0], sol.y[1], self.get_acceleration(sol.y[0])]


    # Solve equations of state explicitly
    # y0: initial displacement [rad]
    # return: [time [s], displacement [rad], velocity [rad/s], acceleration [rad/s/s]]
    def explicit(self, y0, t):
        freq = np.sqrt(-Constants.g_acceleration / self.len)
        return [t, y0*np.cos(freq*t), -y0*freq*np.sin(freq*t), -y0*freq*freq*np.cos(freq*t)]


# Double pendulum: Each part is of the same mass and length. Centre of mass is halfway along the length.
class DoublePendulum:
    # Constructor
    # length: length [m]
    # mass:   mass of each rod [kg]
    def __init__(self, length = 1.0, mass = 1.0):
        self.len  = length
        self.mass = mass


    # Get the angular velocity of the first pendulum
    # theta1: theta 1 [rad]
    # theta2: theta 2 [rad]
    # p1: dL/d(theta1/dt)
    # p2: dL/d(theta2/dt)
    # return derivative of theta 1
    def __get_theta1_dot(self, theta1, theta2, p1, p2):
        a = np.cos(theta1) * np.cos(theta2) + np.sin(theta1) * np.sin(theta2) # cos(theta1-theta2)
        theta1_dot = (2.0*p1 - 3.0*a*p2)/(16.0 - 9.0*a*a)
        theta1_dot *= 6.0/(self.len*self.len*self.mass)
        return theta1_dot


    # Get the angular velocity of the second pendulum
    # theta1: theta 1 [rad]
    # theta2: theta 2 [rad]
    # p1: dL/d(theta1/dt)
    # p2: dL/d(theta2/dt)
    # return derivative of theta 2
    def __get_theta2_dot(self, theta1, theta2, p1, p2):
        a = np.cos(theta1) * np.cos(theta2) + np.sin(theta1) * np.sin(theta2) # cos(theta1-theta2)
        theta2_dot =  (8.0*p2 - 3.0*a*p1)/(16.0 - 9.0*a*a)
        theta2_dot *= 6.0/(self.len*self.len*self.mass)
        return theta2_dot


    # Get the kinetic energy of the system
    # theta1: theta 1 [rad]
    # theta2: theta 2 [rad]
    # theta1_dot: d(theta1)/dt
    # theta2_dot: d(theta2)/dt
    # return the kinetic energy of the system in joules
    def get_kinetic_energy(self, theta1, theta2, theta1_dot, theta2_dot):
        a = np.cos(theta1)*np.cos(theta2) + np.sin(theta1)*np.sin(theta2)   # cos(theta1-theta2)
        return (self.mass*self.len*self.len/6.0) * (theta1_dot*theta1_dot + 4.0*theta1_dot*theta1_dot + 3.0*theta1_dot*theta2_dot*a)


    # Get the potential energy of the system
    # theta1: theta 1 [rad]
    # theta2: theta 2 [rad]
    # theta1_dot: d(theta1)/dt
    # theta2_dot: d(theta2)/dt
    # return the potential energy of the system in joules
    def get_potential_energy(self, theta1, theta2, theta1_dot, theta2_dot):
        return -0.5*self.mass*Constants.g_acceleration*(3.0*np.cos(theta1) + np.cos(theta2))


    # Get the total energy of the system
    # theta1: theta 1 [rad]
    # theta2: theta 2 [rad]
    # theta1_dot: d(theta1)/dt
    # theta2_dot: d(theta2)/dt
    # return the total energy of the system in joules
    def get_total_energy(self, theta1, theta2, theta1_dot, theta2_dot):
        return self.get_kinetic_energy  (theta1, theta2, theta1_dot, theta2_dot) \
             + self.get_potential_energy(theta1, theta2, theta1_dot, theta2_dot)


    # Equations of state for the double pendulum
    # t: time [seconds]
    # x: [theta 1 [rad], theta 2 [rad], dL/d(dtheta1/dt)) [second/rad], dL/d(dtheta2/dt)) [second, rad]]
    # return derivative of x
    def deriv(self, t, x):
        theta1 = x[0]   # First  angle in radians
        theta2 = x[1]   # Second angle in radians
        p1     = x[2]   # Partial derivative of the lagrangian with respect to the time derivative of the first  angle
        p2     = x[3]   # Partial derivative of the lagrangian with respect to the time derivative of the second angle

        theta1_dot = self.__get_theta1_dot(theta1, theta2, p1, p2)
        theta2_dot = self.__get_theta2_dot(theta1, theta2, p1, p2)

        dL_dtheta1 =  theta1_dot*theta2_dot*np.sin(theta1-theta2) + (-3.0*Constants.g_acceleration/self.len) * np.sin(theta1)
        dL_dtheta2 = -theta1_dot*theta2_dot*np.sin(theta1-theta2) + (-Constants.g_acceleration/self.len) * np.sin(theta2)

        dL_dtheta1 *= -0.5*self.mass*self.len*self.len
        dL_dtheta2 *= -0.5*self.mass*self.len*self.len

        return [theta1_dot, theta2_dot, dL_dtheta1, dL_dtheta2]


    # Solve equations of state numerically
    # y0: initial state [theta 1 [rad], theta 2 [rad], dL/d(dtheta1/dt)) [second/rad], dL/d(dtheta2/dt)) [second, rad]]
    # t:  time [s]
    # return: [time [s], total energy [J], kinetic energy [J], potential energy [J],
    #          theta 1 [rad], theta 2 [rad], d(theta1)/dt [rad/s], d(theta2)/dt [rad/s]]
    def solve(self, y0, t):
        sol = ivp(self.deriv, [t[0], t[-1]], y0, t_eval=t)

        # Constrain angles to range [-pi,+pi]
        theta1 = sol.y[0] % (2.0 * np.pi)
        theta2 = sol.y[1] % (2.0 * np.pi)

        theta1 -= 2.0*np.pi*(theta1 > np.pi)
        theta2 -= 2.0*np.pi*(theta2 > np.pi)

        theta1_dot = self.__get_theta1_dot(theta1, theta2, sol.y[2], sol.y[3])
        theta2_dot = self.__get_theta2_dot(theta1, theta2, sol.y[2], sol.y[3])

        T = self.get_kinetic_energy  (theta1, theta2, theta1_dot, theta2_dot)
        V = self.get_potential_energy(theta1, theta2, theta1_dot, theta2_dot)
        E = T + V

        return [sol.t, E, T, V, theta1, theta2, theta1_dot, theta2_dot]
