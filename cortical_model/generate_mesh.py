#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
Cortical Bone – Geometry and Mesh Generation
===============================================================================

This script defines the geometric parameters of the cortical bone model,
generates the parametric geometry, calculates the porosity, and generates
the corresponding Gmsh mesh.

The geometry generation functions are defined in cortical.py.

Authors
-------
María José Cervantes
Ramiro M. Irastorza
===============================================================================
"""

from cortical_geometry import *


# =============================================================================
# USER INPUT
# =============================================================================

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    # Sample
    # -------------------------------------------------------------------------

    sample_name = "Poro10_5"

    SAMPLE_parameters.cube_size = 1.0      # [mm]
    SAMPLE_parameters.density = 10         # [1/mm²]


    # -------------------------------------------------------------------------
    # Osteons
    # -------------------------------------------------------------------------

    OSTEON_parameters.d_min = 0.10         # [mm]
    OSTEON_parameters.d_max = 0.25         # [mm]


    # -------------------------------------------------------------------------
    # Haversian canals
    # -------------------------------------------------------------------------

    HAVERSIAN_CANALS_parameters.d_min = 0.04   # [mm]
    HAVERSIAN_CANALS_parameters.d_max = 0.09   # [mm]

    HAVERSIAN_CANALS_parameters.a_min = 0      # [degree]
    HAVERSIAN_CANALS_parameters.a_max = 15     # [degree]


    # -------------------------------------------------------------------------
    # Volkmann canals
    # -------------------------------------------------------------------------

    VOLKMANNS_CANALS_parameters.d_min = 0.04   # [mm]
    VOLKMANNS_CANALS_parameters.d_max = 0.06   # [mm]

    VOLKMANNS_CANALS_parameters.a_min = 0      # [degree]
    VOLKMANNS_CANALS_parameters.a_max = 15     # [degree]


    # -------------------------------------------------------------------------
    # Random seed
    # -------------------------------------------------------------------------

    Seed.semilla = 5


    # =========================================================================
    # GENERATE GEOMETRY
    # =========================================================================

    geometry = generate_model(
        sample_name,
        ".",
    )


    # =========================================================================
    # CALCULATE POROSITY
    # =========================================================================

    porosity = calculate_porosity(
        sample_name,
        ".",
    )


    # =========================================================================
    # GENERATE MESH
    # =========================================================================

    generate_mesh(
        sample_name,
        ".",
    )


    # =========================================================================
    # RESULTS
    # =========================================================================

    print()
    print("---------------------------------------------")
    print("MODEL GENERATION COMPLETED")
    print("---------------------------------------------")
    print(f"Sample             : {sample_name}")
    print(f"Haversian canals   : {geometry[1]}")
    print(f"Volkmann canals    : {geometry[0]}")
    print(f"Porosity           : {100 * porosity:.2f} %")
    print(f"Geometry           : {sample_name}.geo")
    print(f"Mesh               : {sample_name}.msh")