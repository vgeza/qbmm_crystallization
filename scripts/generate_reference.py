#!/usr/bin/env python

import numpy as np
from scipy import interpolate
from scipy.stats import lognorm
from os.path import join

dirname = "reference"

L0 = 2e-8
L1 = 1e-6
L2 = 200e-6
N = 199 * 2**3 + 1
f_max = 8.672e13

mu = 100e-6
sigma = 13.5e-6
# Convert from mean and standard deviation to lognormal distribution parameters
# See e.g. https://se.mathworks.com/help/stats/lognormal-distribution.html
sigma_log = np.sqrt(np.log(1 + (sigma / mu) ** 2))
mu_log = 2 * np.log(mu) - 0.5 * np.log(sigma**2 + mu**2)
fun = lognorm(s=sigma_log, scale=np.exp(mu_log))
print(f"{fun.mean()=:e} {fun.std()=:e}")

L = np.linspace(L1, L2, N)
f = fun.pdf(L)
f *= f_max / f.max()

T_data = np.loadtxt("constant/T-curve.dat")
T_curve = interpolate.interp1d(T_data[:, 0], T_data[:, 1])


def calculate_moment(L, f, order=0):
    m = 0
    dL = np.diff(L)
    LMid = (L[1:] + L[:-1]) / 2
    fMid = (f[1:] + f[:-1]) / 2
    m = np.sum(fMid * LMid**order * dL)
    return m


def calculate_T(t):
    return T_curve(t)


def calculate_G(T, L, C):
    C_sat = calculate_C_sat(T)
    dC = np.clip(C - C_sat)
    return (
        kG
        * np.exp(-Eg / (R * T))
        * np.pow(1.0 + alphaG * L, Lg)
        * np.pow(dC, g)
    )


def calculate_B(T, L, C, m3):
    C_sat = calculate_C_sat(T)
    dC = np.clip(C - C_sat)
    M = kv * rho * m3
    return kS * np.pow(dC, s) * np.exp(-Es / (R * T)) * np.pow(M, Ms)


def calculate_dG_dL(T, L, C):
    if Lg <= 0:
        return 0
    C_sat = calculate_C_sat(T)
    dC = np.clip(C - C_sat)
    return (
        kG
        * np.exp(-Eg / (R * T))
        * np.pow(1.0 + alphaG * L, Lg - 1)
        * alphaG
        * np.pow(dC, g)
    )


def calculate_C_sat(T):
    TC = T - 273.15
    return (
        4e-5 * np.pow(TC, 4)
        - 0.0034 * np.pow(TC, 3)
        + 0.1024 * np.pow(TC, 2)
        - 0.6255 * TC
        + 13.237
    )


def read_param(name, default_value=0):
    files = [
        "constant/crystallizationProperties",
        "constant/populationBalanceProperties",
    ]
    for file in files:
        with open(file) as f:
            for line in f:
                line = line.strip()
                if line.startswith(f"{name} "):
                    # print(line)
                    return float(line.split()[-1][:-1])
    return default_value


def write_moment(order, value):
    filename = join("0", f"moment.{order}.populationBalance")
    lines = open(filename).readlines()

    with open(filename, "w") as f_out:
        for line in lines:
            if line.startswith("internalField"):
                line = f"internalField   uniform {value:.8e};\n"
            f_out.write(line)


kG = read_param("kG", 1e-8)
g = read_param("g", 0)
Eg = read_param("Eg", 0)
Lg = read_param("Lg", 0)
alphaG = read_param("alphaG", 0)
R = read_param("R", 8.3145)
rho = read_param("rhoCrystal", 1200)
kv = read_param("kv", np.pi / 6)

kS = read_param("kS", 0)
s = read_param("s", 1)
Ms = read_param("Ms", 2)
Es = read_param("Es", 0)

C = 40.0
t = 0
dt = 4
t_max = 9000
orders = range(6)

with open(join(dirname, "moments.dat"), "w") as file_out:
    file_out.write("t\tT\tC_sat\tC\tB")
    for order in orders:
        file_out.write(f"\tM{order}")
    file_out.write("\n")
    # output results at t=0
    T = calculate_T(t)
    C_sat = calculate_C_sat(T)
    m3old = calculate_moment(L, f, 3)
    B = calculate_B(T, L, C, m3old)
    # first, output results at zero time
    file_out.write(f"{t}\t{T}\t{C_sat}\t{C}\t{B}")
    np.savetxt(join(dirname, f"CSD-t{t}.dat"), np.column_stack((L, f)))
    for order in orders:
        m = calculate_moment(L, f, order)
        file_out.write(f"\t{m:.9e}")
        write_moment(order, m)
    file_out.write("\n")
    while t < t_max:
        t += dt
        T = calculate_T(t)
        C_sat = calculate_C_sat(T)
        G = calculate_G(T, L, C)
        dGdL = calculate_dG_dL(T, L, C)
        B = calculate_B(T, L, C, m3old)
        # method of characteristics
        L += G * dt
        f -= f * dGdL * dt
        if L[0] > L0:
            L = np.append([0], L)
            f = np.append([0], f)
        f[0] = B / calculate_G(T, 0, C)
        if int(t) % 300 == 0:
            np.savetxt(join(dirname, f"CSD-t{t}.dat"), np.column_stack((L, f)))
        m3 = calculate_moment(L, f, 3)
        C -= kv * rho * (m3 - m3old)
        file_out.write(f"{t}\t{T}\t{C_sat}\t{C}\t{B}")
        for order in orders:
            m = calculate_moment(L, f, order)
            if order == 0:
                negative_m0 = m < 0
            file_out.write(f"\t{m:.9e}")
        file_out.write("\n")
        m3old = m3
        if negative_m0:
            print("M0 is negative, stopping simulation")
            break
