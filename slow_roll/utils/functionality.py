import numpy as np

from scipy.integrate import odeint
import matplotlib.pyplot as plt

# from utils.constants import (
#     VO,
#     TIME_UNIT,
#     a_0,
# )

V0 = 22
TIME_UNIT = 5e-5
a_0 = 1e-3


def set_hubble_parameter(
    phi_0,
    phi_prime,
):
    return np.sqrt(phi_prime**2 / 6 + (V0 * (phi_0**2) / (3 * TIME_UNIT**2)))


def set_values():
    phi_0 = 1
    phi_prime_0 = 0
    h_0 = set_hubble_parameter(phi_0, phi_prime_0)
    print("set_values: ", [phi_0, phi_prime_0, h_0, a_0])
    return [phi_0, phi_prime_0, h_0, a_0]


def set_system_diff_eqs(var, T):
    [x, y, z, A] = var

    # Note that all derivatives are taken wrt the scaled, dimenstionless cosmic time T

    dxdT = y
    dydT = -3 * z * y - V0 * (2 * x) / TIME_UNIT**2
    dzdT = -0.5 * y**2  # = -z**2 + (v0*f(x)/S**2 - y**2)/3
    dAdT = A * z
    print("set_system_diff_eqs: ", [dxdT, dydT, dzdT, dAdT])
    return [dxdT, dydT, dzdT, dAdT]


def solve_klein(T):
    sol = odeint(set_system_diff_eqs, set_values(), T, mxstep=90)

    x, y, z, A = np.transpose(sol)
    phi, phi_t, H = x, y * TIME_UNIT, z * TIME_UNIT
    print("LLEGA HASTA ACÁ")
    return phi, phi_t, H


def plot_solution():
    T = np.linspace(0, 1, 10)
    phi, phi_t, H = solve_klein(T)
    plt.plot(phi, phi_t, "r", lw=2)
    plt.show()


plot_solution()
