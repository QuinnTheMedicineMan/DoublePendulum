#!/usr/bin/env python

import numpy as np
from scipy.integrate import solve_ivp as ivp
import matplotlib.pyplot as plt
import Constants


# Simple Harmonic Pendulum
class SimpleHarmonicPendulum:
    # Constructor
    # length: length of pendulum [m]
    def __init__(self, length = 1.0):
        self.len = length


    # Get simple harmonic oscillator acceleration
    # x: displacement [m]
    # return: acceleration [m/s]
    def get_acceleration(self, x):
        return -self.len * x


    # Call operator: Equation of state for the SHO
    # t: time [seconds]
    # x: [displacement [m], velocity [m/s]]
    # return: [velocity [m/s], acceleration [m/s/s]]
    def __call__(self, t, x):
        return [x[1], self.get_acceleration(x[0])]


# Double pendulum: Each part is of the same mass and length. Centre of mass is halfway along the length.
class DoublePendulum:
    # Constructor
    # len: length [m]
    def __init__(self, len = 1.0):
        self.len = len

    # Call operator: Equations of state for the double pendulum
    def __call__(self, t, x):
