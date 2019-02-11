#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import SimpleHarmonicOscillator as sho


def main():

    # Solve numerically
    t   = np.linspace(0, 30, 3000)
    pen = sho.SimpleHarmonicPendulum()

    sol = pen.solve([1.0, 0.0], t)

    plt.plot(sol[0], np.degrees(sol[1]), label="Displacement [deg]")
    plt.plot(sol[0], np.degrees(sol[2]), label="Velocity     [deg/s]")
    plt.plot(sol[0], np.degrees(sol[3]), label="Acceleration [deg/s/s]")

    plt.title("Simple Harmonic Pendulum")
    plt.xlabel("Time [s]")
    plt.legend()
    plt.show()

    # Calculate explicitly
    exp = pen.explicit([1.0, 0.0], t)

    plt.plot(sol[0], sol[1] - exp[1], label="Displacement error [deg]")
    plt.plot(sol[0], sol[2] - exp[2], label="Velocity error     [deg/s]")
    plt.plot(sol[0], sol[3] - exp[3], label="Acceleration error [deg/s/s]")

    plt.title("Simple Harmonic Pendulum Errors")
    plt.xlabel("Time [s]")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
