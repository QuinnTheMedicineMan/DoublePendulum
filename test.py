#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import Pendulums as sho


def main():


    # Simple harmonic pendulum checks
    # Solve numerically
    t   = np.linspace(0, 30, 3000)
    pen = pn.SimpleHarmonicPendulum()

    sol = pen.solve([np.radians(10.0), 0.0], t)

    plt.plot(sol[0], np.degrees(sol[1]), label="Displacement [deg]")
    plt.plot(sol[0], np.degrees(sol[2]), label="Velocity     [deg/s]")
    plt.plot(sol[0], np.degrees(sol[3]), label="Acceleration [deg/s/s]")

    plt.title("Simple Harmonic Pendulum")
    plt.xlabel("Time [s]")
    plt.legend()
    plt.show()

    # Calculate explicitly
    exp = pen.explicit(np.radians(10.0), t)

    plt.plot(sol[0], sol[1] - exp[1], label="Displacement error [deg]")
    plt.plot(sol[0], sol[2] - exp[2], label="Velocity error     [deg/s]")
    plt.plot(sol[0], sol[3] - exp[3], label="Acceleration error [deg/s/s]")

    plt.title("Simple Harmonic Pendulum Errors")
    plt.xlabel("Time [s]")
    plt.legend()
    plt.show()


    # Double pendulum checks
    # Does a pendulum at rest stay put?
    pen = pn.DoublePendulum()

    sol = pen.solve([0.0, 0.0, 0.0, 0.0], t)

    plt.plot(sol[0], np.degrees(sol[4]), label=r"$\theta_1$ [deg]")
    plt.plot(sol[0], np.degrees(sol[5]), label=r"$\theta_2$ [deg]")

    plt.title("Double Pendulum at Rest")
    plt.xlabel("Time [s]")
    plt.legend()
    plt.show()


    # Make sure total energy of the system is conserved
    plt.plot(sol[0], sol[1], label="Total Energy",     color='black')
    plt.plot(sol[0], sol[2], label="Kinetic Energy",   color='red')
    plt.plot(sol[0], sol[3], label="Potential Energy", color='blue')

    plt.title("Double Pendulum at Rest")
    plt.xlabel("Time [s]")
    plt.ylabel("Energy [J]")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
