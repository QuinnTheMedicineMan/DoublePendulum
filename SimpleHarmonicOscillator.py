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
    # return: [time [s], displacement [rad], velocity [rad/s], acceleration [rad/s/s]]
    def solve(self, y0, t):
        sol = ivp(self.deriv, [t[0], t[-1]], y0, t_eval=t)
        return [sol.t, sol.y[0], sol.y[1], self.get_acceleration(sol.y[0])]


    # Solve equations of state explicitly
    # y0: initial displacement [rad]
    # return: [time [s], displacement [rad], velocity [rad/s], acceleration [rad/s/s]]
    def explicit(self, y0, t):
        freq = np.sqrt(-Constants.g_acceleration / self.len)
        return [t, np.cos(freq*t), np.sin(freq*t), -np.cos(freq*t)]


# Double pendulum: Each part is of the same mass and length. Centre of mass is halfway along the length.
class DoublePendulum:
    # Constructor
    # length: length [m]
    # mass:   mass of each rod [kg]
    def __init__(self, length = 1.0, mass = 1.0):
        self.len  = length
        self.mass = mass


    # Equations of state for the double pendulum
    # t: time [seconds]
    # x: [theta 1 [rad], theta 2 [rad], dL/d(dtheta1/dt)) [second/rad], dL/d(dtheta2/dt)) [second, rad]]
    # return derivative of x
    def deriv(self, t, x):
        theta1 = x[0]   # First  angle in radians
        theta2 = x[1]   # Second angle in radians
        p1     = x[2]   # Partial derivative of the lagrangian with respect to the time derivative of the first  angle
        p2     = x[3]   # Partial derivative of the lagrangian with respect to the time derivative of the second angle

        # Potential for speed up / avoidance of floating point errors in the cosine
        theta1_dot = (2.0*p1 - 3.0*np.cos(theta1-theta2)*p2) / (16.0 - 9*np.cos(theta1-theta2)*np.cos(theta1-theta2))
        theta2_dot = (8.0*p2 - 3.0*np.cos(theta1-theta2)*p1) / (16.0 - 9*np.cos(theta1-theta2)*np.cos(theta1-theta2))

        theta1_dot *= 6.0/(self.len*self.len*self.mass)
        theta2_dot *= 6.0/(self.len*self.len*self.mass)

        dL_dtheta1 =  theta1_dot*theta2_dot*np.sin(theta1-theta2) + (3.0*Constants.g_acceleration/self.len) * np.sin(theta1)
        dL_dtheta2 = -theta1_dot*theta2_dot*np.sin(theta1-theta2) + (Constants.g_acceleration/self.len) * np.sin(theta2)

        dL_dtheta1 *= -0.5*self.mass*self.len*self.len
        dL_dtheta2 *= -0.5*self.mass*self.len*self.len

        return [theta1_dot, theta2_dot, dL_dtheta1, dL_dtheta2]


