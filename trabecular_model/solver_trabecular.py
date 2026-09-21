#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
===============================================================================
FEniCSx – Effective Electrical Conductivity of Trabecular Bone
===============================================================================

This script computes the effective electrical conductivity of three-dimensional 
in silico cylindrical trabecular bone samples based on micro-CT reconstructed 
geometries.

The in silico samples are represented as heterogeneous two-phase domains
composed of:

    - Bone matrix  (physical tag 200)
    - Marrow space (physical tag 100)

The electrical potential distribution is obtained by solving the steady-state
conduction equation

        div(sigma * grad(V)) = 0

using the finite element method (FEM) implemented in FEniCSx.

Dirichlet boundary conditions are imposed on the two opposite surfaces of the
cylindrical samples, identified by physical facet tags in the Gmsh meshes.

The effective electrical conductivity is calculated from the total electrical
power dissipated in the domain according to

        R = DeltaV^2 / P

        sigma_eff = L / (R * A)

where L is the sample length, A is the cross-sectional area, DeltaV is the
imposed potential difference, and P is the total dissipated electrical power.

The micro-CT-based geometries are defined in millimeters, while all FEM
calculations are performed using SI units.

Units
-----
Input geometry                  : mm
Internal FEM geometry           : m
Input conductivity              : mS/m
Internal FEM conductivity       : S/m
Effective conductivity output   : mS/m

Requirements
------------
Python      : 3.12.3
DOLFINx     : 0.9.0
petsc4py    : 3.19.6
mpi4py      : 3.1.5
NumPy       : 1.26.4
UFL         : 2024.2.0

Authors
-------
María José Cervantes
Ramiro M. Irastorza

===============================================================================
"""


# =============================================================================
# IMPORTS
# =============================================================================

import os
import numpy as np
import ufl
from petsc4py import PETSc


from mpi4py import MPI
from dolfinx import mesh, fem, default_scalar_type
from dolfinx.io import gmshio
from dolfinx.fem.petsc import LinearProblem


# =============================================================================
# FEM SOLVER
# =============================================================================

def solve_trabecular_bone(
        name,
        sample_length,
        sample_radius,
        sigma_matrix,
        sigma_marrow):
    """
    Compute the effective electrical conductivity of a cylindrical
    trabecular bone sample using the finite element method.

    Parameters
    ----------
    name : str
        Name of the Gmsh mesh file without the ".msh" extension.

    sample_length : float
        Length of the cylindrical sample in mm.

    sample_radius : float
        Radius of the cylindrical sample in mm.

    sigma_matrix : float
        Electrical conductivity of the bone matrix in mS/m.

    sigma_marrow : float
        Electrical conductivity of the marrow phase in mS/m.

    Returns
    -------
    float
        Effective electrical conductivity in mS/m.

    Notes
    -----
    Physical tags used in the mesh:

        200 : bone matrix
        100 : marrow space

    The longitudinal axis of the cylindrical sample is assumed to be
    aligned with the z-axis.
    """

    # -------------------------------------------------------------------------
    # Geometry
    # -------------------------------------------------------------------------

    msh_file = os.path.join(f"{name}.msh")

    # Convert sample dimensions from mm to m
    length = sample_length * 1e-3
    radius = sample_radius * 1e-3

    # Cross-sectional area [m²]
    area = np.pi * radius**2


    # -------------------------------------------------------------------------
    # Electrical properties
    # -------------------------------------------------------------------------

    # Convert conductivities from mS/m to S/m
    sigma_bone = sigma_matrix * 1e-3
    sigma_mar = sigma_marrow * 1e-3


    # -------------------------------------------------------------------------
    # Applied potential
    # -------------------------------------------------------------------------

    U0 = 40.0
    Ug = 0.0

    delta_V = U0 - Ug


    # -------------------------------------------------------------------------
    # Read mesh
    # -------------------------------------------------------------------------

    msh, cell_tags, facet_tags = gmshio.read_from_msh(
        msh_file,
        MPI.COMM_WORLD,
        0,
        gdim=3
    )

    # Convert mesh coordinates from mm to m
    msh.geometry.x[:] *= 1e-3

    tdim = msh.topology.dim
    fdim = tdim - 1

    msh.topology.create_connectivity(fdim, tdim)


    # -------------------------------------------------------------------------
    # Function space
    # -------------------------------------------------------------------------

    V = fem.functionspace(msh, ("P", 1))


    # -------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------
    
    # Electrode surfaces are identified by physical facet tags
    # defined in the Gmsh mesh:
    #
    #     10 : active electrode
    #     20 : passive electrode
    
    active_facets = facet_tags.find(10)
    passive_facets = facet_tags.find(20)
    
    print("Active electrode facets (tag 10):", len(active_facets))
    print("Passive electrode facets (tag 20):", len(passive_facets))
    
    
    # Locate degrees of freedom on the electrode surfaces
    
    active_dofs = fem.locate_dofs_topological(
        V,
        fdim,
        active_facets,
    )
    
    passive_dofs = fem.locate_dofs_topological(
        V,
        fdim,
        passive_facets,
    )
    
    
    # Apply Dirichlet boundary conditions
    
    bc_active = fem.dirichletbc(
        PETSc.ScalarType(U0),
        active_dofs,
        V,
    )
    
    bc_passive = fem.dirichletbc(
        PETSc.ScalarType(Ug),
        passive_dofs,
        V,
    )
    # -------------------------------------------------------------------------
    # Spatial conductivity distribution
    # -------------------------------------------------------------------------

    Q = fem.functionspace(msh, ("DG", 0))

    sigma = fem.Function(Q)

    # Physical tags defined in the Gmsh mesh
    MARROW_TAG = 100
    BONE_MATRIX_TAG = 200
    
  
    marrow_cells = cell_tags.find(MARROW_TAG)
    matrix_cells = cell_tags.find(BONE_MATRIX_TAG)
    
    sigma.x.array[marrow_cells] = sigma_mar
    sigma.x.array[matrix_cells] = sigma_bone

    # -------------------------------------------------------------------------
    # Variational problem
    # -------------------------------------------------------------------------

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    a = (
        ufl.inner(
            sigma * ufl.grad(u),
            ufl.grad(v)
        )
        * ufl.dx
    )

    L = (
        fem.Constant(
            msh,
            default_scalar_type(0)
        )
        * v
        * ufl.dx
    )


    # -------------------------------------------------------------------------
    # Linear solver
    # -------------------------------------------------------------------------

    problem = LinearProblem(
        a,
        L,
        bcs=[bc_active, bc_passive],
        petsc_options={
            "ksp_type": "cg",
            "pc_type": "hypre",
            "ksp_rtol": 1e-10,
            "ksp_atol": 1e-15,
        },
    )

    potential = problem.solve()


    # -------------------------------------------------------------------------
    # Electrical power
    # -------------------------------------------------------------------------

    power = fem.assemble_scalar(
        fem.form(
            sigma
            * ufl.inner(
                ufl.grad(potential),
                ufl.grad(potential)
            )
            * ufl.dx
        )
    )


    # -------------------------------------------------------------------------
    # Effective electrical conductivity
    # -------------------------------------------------------------------------

    resistance = delta_V**2 / power

    sigma_eff = length / (resistance * area)

    # Convert S/m to mS/m
    sigma_eff_mS_m = sigma_eff * 1e3


    return sigma_eff_mS_m


# =============================================================================
# USER INPUT
# =============================================================================
#
# Modify ONLY the parameters in this section.
#
# =============================================================================

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    # Mesh
    # -------------------------------------------------------------------------

    # Mesh file name without the ".msh" extension
    sample_name = "modelo_muestra_3c_231"


    # -------------------------------------------------------------------------
    # Sample dimensions [mm]
    # -------------------------------------------------------------------------

    sample_length = 2.6
    sample_radius = 2.1


    # -------------------------------------------------------------------------
    # Material conductivities [mS/m]
    # -------------------------------------------------------------------------

    sigma_matrix = 20.0
    sigma_marrow = 301.0


    # =========================================================================
    # RUN FEM MODEL
    # =========================================================================

    sigma_eff = solve_trabecular_bone(
        name=sample_name,
        sample_length=sample_length,
        sample_radius=sample_radius,
        sigma_matrix=sigma_matrix,
        sigma_marrow=sigma_marrow
    )


    # =========================================================================
    # RESULTS
    # =========================================================================

    print("\n====================================================")
    print("TRABECULAR BONE FEM RESULTS")
    print("====================================================")

    print(f"Sample                    : {sample_name}")
    print(f"Sample length             : {sample_length:.3f} mm")
    print(f"Sample radius             : {sample_radius:.3f} mm")

    print("----------------------------------------------------")

    print(f"Bone matrix conductivity  : {sigma_matrix:.3f} mS/m")
    print(f"Marrow conductivity       : {sigma_marrow:.3f} mS/m")

    print("----------------------------------------------------")

    print(f"Effective conductivity    : {sigma_eff:.3f} mS/m")

    print("====================================================\n")
