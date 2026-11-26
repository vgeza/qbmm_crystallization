#!/usr/bin/env python

import argparse
import numpy as np
import pint
import yaml
from scipy import interpolate
from scipy.stats import lognorm
from os import makedirs
from os.path import join

VERSION = "1.11"
UREG = pint.UnitRegistry(autoconvert_offset_to_baseunit=True)
Q_ = UREG.Quantity

parser = argparse.ArgumentParser(
    prog="generate_reference",
    description=(
        "Generates an initial state and a reference solution for PBM."
    ),
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-v", "--version", action="version", version=f"%(prog)s {VERSION}"
)

parser.add_argument(
    "-s", "--source", default="parameters.yaml", help="Input YAML file"
)

args = parser.parse_args()

try:
    with open(args.source) as f:
        print(f"Reading {args.source}")
        parameters = yaml.safe_load(f)
except:
    parameters = {}
if not parameters:
    parameters = {}
print(f"\n{parameters=}")
parameters_CSD = parameters.get("CSD", {})
parameters_conditions = parameters.get("conditions", {})
parameters_other = parameters.get("other", {})


def get_T_curve_data(parameters):
    par = parameters.get("conditions", {})
    T0 = Q_(par.get("T0", "20 degC")).to("K")
    t = [0]
    T = [T0.magnitude]
    if "T1" in par:
        T1 = Q_(par.get("T1")).to("K")
        beta0 = Q_(par.get("beta0")).to("K/s")
        t1 = (T1 - T0) / beta0
        T.append(T1.magnitude)
        t.append(t1.magnitude)
    if "T2" in par:
        T2 = Q_(par.get("T2")).to("K")
        beta1 = Q_(par.get("beta1")).to("K/s")
        t2 = t1 + (T2 - T1) / beta1
        T.append(T2.magnitude)
        t.append(t2.magnitude)
    return t, T


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
    dC = np.maximum(C - C_sat, 0)
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
        try:
            with open(file) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(f"{name} "):
                        return float(line.split()[-1][:-1])
        except:
            continue
    return default_value


def write_T_curve(t, T, fname):
    with open(fname, "w") as f:
        for ti, Ti in zip(t, T):
            f.write(f"{ti:g} {Ti:g}\n")


def write_moment(order, value):
    filename = join("0", f"moment.{order}.populationBalance")
    write_value(filename, f"{value:.9e}")


def write_value(filename, value, keyword="internalField"):
    lines = open(filename).readlines()

    with open(filename, "w") as f_out:
        for line in lines:
            if line.startswith(keyword):
                # only change the last bit of the line
                i = line.rfind(" ")
                line = line[:i] + f" {value};\n"
            f_out.write(line)


dirname = parameters_other.get("output", "reference")
makedirs(dirname, exist_ok=True)

t, T = get_T_curve_data(parameters)
T_curve = interpolate.interp1d(
    t, T, bounds_error=False, fill_value=(T[0], T[-1])
)
fname = parameters_other.get("T_curve", "constant/T-curve.dat")
write_T_curve(t, T, fname)

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

C = Q_(parameters_conditions.get("C0", "40 kg/m^3")).to("kg/m^3").magnitude
write_value(join("0", "C"), C)
t = 0
dt = Q_(parameters_other.get("dt", "4 s")).to("s").magnitude
t_max = Q_(parameters_conditions.get("t_max", "150 min")).to("s").magnitude
# keep OpenFOAM time step independent, and write only end time
write_value(join("system", "controlDict"), t_max, "endTime")

orders = range(6)

# generate the initial CSD
mu = Q_(parameters_CSD.get("mu", "100 um")).to("m").magnitude
sigma = Q_(parameters_CSD.get("sigma", "13.5 um")).to("m").magnitude
# Convert from mean and standard deviation to lognormal distribution parameters
# See e.g. https://se.mathworks.com/help/stats/lognormal-distribution.html
sigma_log = np.sqrt(np.log(1 + (sigma / mu) ** 2))
mu_log = 2 * np.log(mu) - 0.5 * np.log(sigma**2 + mu**2)
fun = lognorm(s=sigma_log, scale=np.exp(mu_log))
print(f"{fun.mean()=:e} {fun.std()=:e}")

L0 = Q_(parameters_CSD.get("L0", "0.02 um")).to("m").magnitude
L1 = Q_(parameters_CSD.get("L1", "1 um")).to("m").magnitude
L2 = Q_(parameters_CSD.get("L2", "200 um")).to("m").magnitude
N = int(parameters_CSD.get("N", 199 * 2**3 + 1))
L = np.linspace(L1, L2, N)
f = fun.pdf(L)

for pp in [10, 50, 90]:
    print(f"D{pp}={1e6*fun.ppf(pp/100):g} um")

# recalculate seed loading for unit volume and scale the CSD
m_seed = Q_(parameters_conditions.get("m_seed", "1.28 g")).to("kg")
volume = Q_(parameters_conditions.get("volume", "350 mL")).to("m^3")
m_seed_scaled = m_seed / volume * Q_("1 m^3")
m3 = calculate_moment(L, f, 3)  # fun.moment(3)
scale = m_seed_scaled.magnitude / (kv * rho * m3)
print(f"{scale=:e}")
f *= scale

with open(join(dirname, "moments.dat"), "w") as file_out:
    file_out.write("t\tT\tC_sat\tC\tB\tG0\tf0")
    for order in orders:
        file_out.write(f"\tM{order}")
    file_out.write("\n")
    # output results at t=0
    T = calculate_T(t)
    C_sat = calculate_C_sat(T)
    m3old = calculate_moment(L, f, 3)
    B = calculate_B(T, L, C, m3old)
    G0 = calculate_G(T, 0, C)
    # first, output results at zero time
    file_out.write(f"{t}\t{T}\t{C_sat}\t{C}\t{B}\t{G0}\t{f[0]}")
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
        G0 = calculate_G(T, 0, C)
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
            fname = f"CSD-t{t:g}.dat"
            print(f"Writing {fname}")
            np.savetxt(join(dirname, fname), np.column_stack((L, f)))
        m3 = calculate_moment(L, f, 3)
        C -= kv * rho * (m3 - m3old)
        file_out.write(f"{t}\t{T}\t{C_sat}\t{C}\t{B}\t{G0}\t{f[0]}")
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
        if np.any(
            [
                np.any(np.isnan(G)),
                np.any(np.isnan(dGdL)),
                np.any(np.isnan(B)),
                np.any(np.isnan(f[0])),
                np.any(np.isnan(C)),
            ]
        ):
            print("NaN value encountered, stopping simulation")
            break
