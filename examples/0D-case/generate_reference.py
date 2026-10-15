#!/usr/bin/env python

import numpy as np
from scipy import interpolate
from scipy.stats import lognorm
from os.path import join

dirname = "reference"

L1 = 1e-6
L2 = 200e-6
N = 199 * 2**3 + 1
f_max = 4.55161e16

mu = 100e-6
sigma = 5e-6
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
R = read_param("R", 8.314)
rho = read_param("rhoCrystal", 1200)
kv = read_param("kv", np.pi / 6)

C = 0.6 * rho
t = 0
dt = 1
t_max = 2000
orders = range(6)

with open(join(dirname, "moments.dat"), "w") as file_out:
    file_out.write("t\tT\tC_sat\tC")
    for order in orders:
        file_out.write(f"\tM{order}")
    file_out.write("\n")
    # output results at t=0
    T = calculate_T(t)
    C_sat = calculate_C_sat(T)
    # first, output results at zero time
    file_out.write(f"{t}\t{T}\t{C_sat}\t{C}")
    np.savetxt(join(dirname, f"CSD-t{t}.dat"), np.column_stack((L, f)))
    m3old = calculate_moment(L, f, 3)
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
        # method of characteristics
        L += G * dt
        f -= f * dGdL * dt
        if int(t) % 200 == 0:
            np.savetxt(join(dirname, f"CSD-t{t}.dat"), np.column_stack((L, f)))
        m3 = calculate_moment(L, f, 3)
        C -= kv * rho * (m3 - m3old)
        file_out.write(f"{t}\t{T}\t{C_sat}\t{C}")
        for order in orders:
            m = calculate_moment(L, f, order)
            file_out.write(f"\t{m:.9e}")
        file_out.write("\n")
        m3old = m3
