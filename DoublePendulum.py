#!/usr/bin/env python

import sys
import argparse
import Pendulums as pn
import numpy as np
import matplotlib.pyplot as plt

def main(argv):

    zap = argparse.ArgumentParser()
    zap.add_argument('-t', '--theta',   nargs=2, type=int,   required=True, help='Starting angles for the double pendulum in degrees')
    zap.add_argument('-T', '--trange',  nargs=2, type=float, default=[0.0, 60.0], help='Time range in seconds [min, max]')
    zap.add_argument('-N', '--npoints', type=int,   default=6000, help='Number of time points to sample')
    zap.add_argument('-m', '--mass',    type=float, default=1.0,  help='Mass of each pendulum in kg')
    zap.add_argument('-L', '--len',     type=float, default=1.0,  help='Length of each pendulum in metres')

    args = zap.parse_args()

    theta1 = np.radians(args.theta[0])
    theta2 = np.radians(args.theta[1])
    length = args.len

    pen  = pn.DoublePendulum(length, args.mass)                         # Initialise the double pendulum
    time = np.linspace(args.trange[0], args.trange[1], args.npoints)    # Time range to run the simulation in

    sol  = pen.solve([theta1, theta2, 0.0, 0.0], time)                  # Run the double pendulum simulation

    # Extract output
    total_energy     = sol[1]
    kinetic_energy   = sol[2]
    potential_energy = sol[3]
    theta1           = np.degrees(sol[4])
    theta2           = np.degrees(sol[5])

    # Plot the angle of each pendulum against time
    plt.plot(time, theta1, label=r"$\theta_1$ [deg]")
    plt.plot(time, theta2, label=r"$\theta_2$ [deg]")

    plt.legend()
    plt.title("Double Pendulum")
    plt.xlabel("Time [s]")
    plt.legend()
    plt.show()

    # Plot one angle against the other
    plt.plot(theta1, theta2)

    plt.title("Double Pendulum")
    plt.xlabel(r"$\theta_1$ [deg]")
    plt.ylabel(r"$\theta_2$ [deg]")
    plt.show()

    # Track the x,y coordinates of the first pendulum
    x =  length * np.sin(np.radians(theta1))
    y = -length * np.cos(np.radians(theta1))

    plt.plot(x, y)

    plt.title("Double Pendulum; First Pendulum Position")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.show()

    # Track the x,y coordinates of the second pendulum
    x += length * np.sin(np.radians(theta2))
    y -= length * np.cos(np.radians(theta2))

    plt.plot(x, y)

    plt.title("Double Pendulum; Second Pendulum Position")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.show()

    # Plot the energy of the system against time (should be constant except for floating point errors etc)
    plt.plot(time, total_energy,     label="Total Energy",     color='black')
    plt.plot(time, kinetic_energy,   label="Kinetic Energy",   color='red')
    plt.plot(time, potential_energy, label="Potential Energy", color='blue')

    plt.title("Double Pendulum")
    plt.xlabel("Time [s]")
    plt.ylabel("Energy [J]")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main(sys.argv)
