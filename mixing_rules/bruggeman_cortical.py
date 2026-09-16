# -*- coding: utf-8 -*-
"""
Bruggeman effective-medium model for cortical bone.

The model estimates the effective electrical conductivity
in the x, y, and z directions for a given BV/TV.

@author: majoc
"""

import numpy as np
from scipy.optimize import root
from scipy.constants import epsilon_0


# =========================================================
# GENERAL PARAMETERS
# =========================================================

freq = 100e3  # Hz
omega = 2 * np.pi * freq


# =========================================================
# ELECTRICAL PROPERTIES
# =========================================================

# Cortical bone matrix
cond_1 = 3.84e-3  # S/m
epsr_1 = 2.28e2

# Bone marrow / vascular channels
cond_2 = 267e-3   # S/m
epsr_2 = 1.11e2


# Complex relative permittivities
epsc_1 = epsr_1 - 1j * cond_1 / (epsilon_0 * omega)

epsc_2 = epsr_2 - 1j * cond_2 / (epsilon_0 * omega)
epsc_3 = epsr_2 - 1j * cond_2 / (epsilon_0 * omega)
epsc_4 = epsr_2 - 1j * cond_2 / (epsilon_0 * omega)


# =========================================================
# DEPOLARIZATION TENSORS
# =========================================================

# Cortical bone matrix
N1 = np.diag([1/3, 1/3, 1/3])

# Vertical channels
N2 = np.diag([1/2, 1/2, 0])

# Horizontal channels along x
N3 = np.diag([0, 1/2, 1/2])

# Horizontal channels along y
N4 = np.diag([1/2, 0, 1/2])


# =========================================================
# CHANNEL DISTRIBUTION
# =========================================================

# Total and vertical vascular porosity (%)
phi = 7.60
phi_vertical = 4.98

# Fraction of vascular porosity associated with vertical channels
p_v = phi_vertical / phi


# =========================================================
# BRUGGEMAN MODEL - FOUR PHASES
# =========================================================

def bruggeman(parameters, *data):

    f1, f2, f3, f4, eps1, eps2, eps3, eps4, N1, N2, N3, N4 = data

    ex = parameters[0] + 1j * parameters[1]
    ey = parameters[2] + 1j * parameters[3]
    ez = parameters[4] + 1j * parameters[5]

    eps_eff = np.diag([ex, ey, ez]).astype(complex)

    I = np.eye(3, dtype=complex)

    eps1_tensor = eps1 * I
    eps2_tensor = eps2 * I
    eps3_tensor = eps3 * I
    eps4_tensor = eps4 * I

    A1 = eps_eff + N1 @ (eps1_tensor - eps_eff)
    A2 = eps_eff + N2 @ (eps2_tensor - eps_eff)
    A3 = eps_eff + N3 @ (eps3_tensor - eps_eff)
    A4 = eps_eff + N4 @ (eps4_tensor - eps_eff)

    T1 = f1 * (eps1_tensor - eps_eff) @ np.linalg.inv(A1)
    T2 = f2 * (eps2_tensor - eps_eff) @ np.linalg.inv(A2)
    T3 = f3 * (eps3_tensor - eps_eff) @ np.linalg.inv(A3)
    T4 = f4 * (eps4_tensor - eps_eff) @ np.linalg.inv(A4)

    residual = T1 + T2 + T3 + T4

    d = np.diag(residual)

    return np.array([
        np.real(d[0]), np.imag(d[0]),
        np.real(d[1]), np.imag(d[1]),
        np.real(d[2]), np.imag(d[2])
    ])


# =========================================================
# EFFECTIVE CONDUCTIVITY
# =========================================================

def calculate_effective_conductivity(BVTV):

    # Volume fractions
    f1 = BVTV
    f2 = (1 - f1) * p_v
    f3 = ((1 - f1) * (1 - p_v)) / 2
    f4 = ((1 - f1) * (1 - p_v)) / 2

    # Initial guess
    x0 = np.array([200, 0, 200, 0, 200, 0], dtype=float)

    # Solve Bruggeman equation
    sol = root(
        bruggeman,
        x0,
        args=(
            f1, f2, f3, f4,
            epsc_1, epsc_2, epsc_3, epsc_4,
            N1, N2, N3, N4
        )
    )

    # Effective complex permittivities
    ex = sol.x[0] + 1j * sol.x[1]
    ey = sol.x[2] + 1j * sol.x[3]
    ez = sol.x[4] + 1j * sol.x[5]

    # Effective conductivities (mS/m)
    sigma_x = -np.imag(ex) * epsilon_0 * omega / 1e-3
    sigma_y = -np.imag(ey) * epsilon_0 * omega / 1e-3
    sigma_z = -np.imag(ez) * epsilon_0 * omega / 1e-3

    return sigma_x, sigma_y, sigma_z


# =========================================================
# EXAMPLE
# =========================================================

# Set the desired bone volume fraction (BV/TV)
BVTV = 0.90

sigma_x, sigma_y, sigma_z = calculate_effective_conductivity(BVTV)

print(f"BV/TV = {BVTV:.3f}")
print(f"Effective conductivity x = {sigma_x:.3f} mS/m")
print(f"Effective conductivity y = {sigma_y:.3f} mS/m")
print(f"Effective conductivity z = {sigma_z:.3f} mS/m")