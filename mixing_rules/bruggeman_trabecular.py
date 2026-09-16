# -*- coding: utf-8 -*-
"""
Bruggeman effective-medium model for trabecular bone.

The model estimates the effective electrical conductivity
for a given BV/TV.

@author: majoc
"""

import numpy as np
from scipy.constants import epsilon_0
from scipy.optimize import newton


# =========================================================
# GENERAL PARAMETERS
# =========================================================

freq = 100e3  # Hz
omega = 2 * np.pi * freq


# =========================================================
# ELECTRICAL PROPERTIES
# =========================================================

# Bone matrix
cond_1 = 2.08e-2  # S/m
epsr_1 = 2.28e2

# Free water
cond_2 = 1.2      # S/m
epsr_2 = 78

# Bone marrow
cond_3 = 1.03e-1  # S/m
epsr_3 = 1.73e2

# Effective relative permittivity
epsr_eff = 4.72e2


# =========================================================
# COMPLEX PERMITTIVITY
# =========================================================

def complex_permittivity(eps_r, conductivity):

    return (
        eps_r
        - 1j * conductivity / (epsilon_0 * omega)
    )


# =========================================================
# BRUGGEMAN MODEL - THREE PHASES
# =========================================================

def bruggeman(cond_eff, f1, f2):

    # Effective medium
    eps_eff = complex_permittivity(
        epsr_eff,
        cond_eff
    )

    # Phase 1: bone matrix
    eps_1 = complex_permittivity(
        epsr_1,
        cond_1
    )

    # Phase 2: free water
    eps_2 = complex_permittivity(
        epsr_2,
        cond_2
    )

    # Phase 3: bone marrow
    eps_3 = complex_permittivity(
        epsr_3,
        cond_3
    )

    # Bone marrow volume fraction
    f3 = 1 - f1 - f2

    residual = (
        f1 * (eps_1 - eps_eff)
        / (eps_1 + 2 * eps_eff)

        + f2 * (eps_2 - eps_eff)
        / (eps_2 + 2 * eps_eff)

        + f3 * (eps_3 - eps_eff)
        / (eps_3 + 2 * eps_eff)
    )

    return np.abs(residual)


# =========================================================
# EFFECTIVE CONDUCTIVITY
# =========================================================

def calculate_effective_conductivity(BVTV):

    # Bone matrix volume fraction
    f1 = BVTV

    # Free-water volume fraction
    f2 = 11 / 30 - BVTV / 3

    # Solve Bruggeman equation
    solution = newton(
        bruggeman,
        x0=2.08e-2,
        args=(f1, f2),
        maxiter=800,
        rtol=0.1
    )

    # Effective conductivity (mS/m)
    sigma = solution / 1e-3

    return sigma


# =========================================================
# EXAMPLE
# =========================================================

# Set the desired bone volume fraction (BV/TV)
BVTV = 0.44

sigma = calculate_effective_conductivity(BVTV)

print(f"BV/TV = {BVTV:.3f}")
print(f"Effective conductivity = {sigma:.3f} mS/m")
