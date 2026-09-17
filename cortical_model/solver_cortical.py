#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
FEniCSx – Effective Electrical Conductivity of Cortical Bone
===============================================================================

This script computes the effective electrical conductivity of three-dimensional
in silico cubic cortical bone samples based on parametrically generated
geometries.

The in silico samples are represented as heterogeneous two-phase domains
composed of:

    - Cortical bone matrix
    - Marrow

The electrical potential distribution is obtained by solving the steady-state
conduction equation

        div(sigma * grad(V)) = 0

using the finite element method (FEM) implemented in FEniCSx.

Dirichlet boundary conditions are imposed on opposite surfaces of the cubic
samples according to their spatial coordinates. The direction of the applied
electric field can be selected along the x-, y-, or z-axis.

The effective electrical conductivity is calculated from the total electrical
power dissipated in the domain according to

        R = DeltaV^2 / P

        sigma_eff = L / (R * A)

where L is the sample length along the applied electric field, A is the
cross-sectional area perpendicular to the field, DeltaV is the imposed potential
difference, and P is the total dissipated electrical power.

The parametrically generated geometries are defined in meters, and all FEM
calculations are performed using SI units.

Units
-----
Input geometry                  : m
Internal FEM geometry           : m
Input conductivity              : mS/m
Internal FEM conductivity       : S/m
Effective conductivity output   : mS/m

Requirements
------------
Python      : 3.10.20
DOLFINx     : 0.7.3
PETSc       : 3.20.6
petsc4py    : 3.20.5
mpi4py      : 3.1.6
NumPy       : 1.26.4
UFL         : 2023.2.0

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

def solve_cortical_bone(
        name,
        sample_length,
        sigma_matrix,
        sigma_marrow,
        direction):
    """
    Compute the effective electrical conductivity of a cubic cortical bone
    sample using the finite element method.

    Parameters
    ----------
    name : str
        Name of the Gmsh mesh file without the ".msh" extension.

    sample_length : float
        Side length of the cubic sample in m.

    sigma_matrix : float
        Electrical conductivity of the cortical bone matrix in mS/m.

    sigma_marrow : float
        Electrical conductivity of the marrow phase in mS/m.

    direction : str
        Direction of the applied electric field: "x", "y", or "z".

    Returns
    -------
    float
        Effective electrical conductivity in mS/m.
    """


    # -------------------------------------------------------------------------
    # Geometry
    # -------------------------------------------------------------------------

    msh_file = os.path.join(f"{name}.msh")

    # Cubic sample
    length = sample_length

    # Cross-sectional area [m²]
    area = sample_length**2


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

    # The cortical mesh is already defined in meters.
    # No geometrical unit conversion is required.

    tdim = msh.topology.dim
    fdim = tdim - 1

    msh.topology.create_connectivity(fdim, tdim)


    # -------------------------------------------------------------------------
    # Function space
    # -------------------------------------------------------------------------

    V = fem.functionspace(msh, ("P", 1))


        # -------------------------------------------------------------------------
    # Selected direction
    # -------------------------------------------------------------------------

    axis_dict = {"x": 0, "y": 1, "z": 2}

    if direction not in axis_dict:
        raise ValueError("direction must be 'x', 'y', or 'z'")

    ax = axis_dict[direction]


    # -------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------

    # Dirichlet boundary conditions are imposed on opposite surfaces
    # of the cubic sample along the selected direction.
    #
    # x[0] = x-axis
    # x[1] = y-axis
    # x[2] = z-axis

    def active_boundary(x):
        return np.isclose(x[ax], sample_length / 2)

    def passive_boundary(x):
        return np.isclose(x[ax], -sample_length / 2)


    active_facets = mesh.locate_entities_boundary(
        msh,
        fdim,
        active_boundary
    )

    passive_facets = mesh.locate_entities_boundary(
        msh,
        fdim,
        passive_boundary
    )


    active_dofs = fem.locate_dofs_topological(
        V,
        fdim,
        active_facets
    )

    passive_dofs = fem.locate_dofs_topological(
        V,
        fdim,
        passive_facets
    )


    bc_active = fem.dirichletbc(
        PETSc.ScalarType(U0),
        active_dofs,
        V
    )

    bc_passive = fem.dirichletbc(
        PETSc.ScalarType(Ug),
        passive_dofs,
        V
    )


    # -------------------------------------------------------------------------
    # Spatial conductivity distribution
    # -------------------------------------------------------------------------

    Q = fem.functionspace(msh, ("DG", 0))

    sigma = fem.Function(Q)


    # Physical tags defined in the Gmsh mesh
    MARROW_TAG = 500
    BONE_MATRIX_TAG = 600
 

    marrow_cells = cell_tags.find(
        MARROW_TAG
    )
    matrix_cells = cell_tags.find(
            BONE_MATRIX_TAG
    )

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
        bcs=[
            bc_active,
            bc_passive
        ],
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

    sigma_eff = length / (
        resistance * area
    )


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
    sample_name = "Poro6"


    # -------------------------------------------------------------------------
    # Direction
    # -------------------------------------------------------------------------

    # Applied electric field direction: "x", "y", or "z"
    direction = "y"


    # -------------------------------------------------------------------------
    # Sample dimensions [m]
    # -------------------------------------------------------------------------

    sample_length = 1.0e-3


    # -------------------------------------------------------------------------
    # Material conductivities [mS/m]
    # -------------------------------------------------------------------------

    sigma_matrix = 3.84
    sigma_marrow = 267.0


    # =========================================================================
    # RUN FEM MODEL
    # =========================================================================

    sigma_eff = solve_cortical_bone(
        name=sample_name,
        sample_length=sample_length,
        sigma_matrix=sigma_matrix,
        sigma_marrow=sigma_marrow,
        direction=direction
    )


    # =========================================================================
    # RESULTS
    # =========================================================================

    print("\n====================================================")
    print("CORTICAL BONE FEM RESULTS")
    print("====================================================")

    print(f"Sample                    : {sample_name}")
    print(f"Direction                 : {direction.upper()}")
    print(f"Sample length             : {sample_length:.6e} m")

    print("----------------------------------------------------")

    print(f"Bone matrix conductivity  : {sigma_matrix:.3f} mS/m")
    print(f"Marrow conductivity       : {sigma_marrow:.3f} mS/m")

    print("----------------------------------------------------")

    print(f"Effective conductivity    : {sigma_eff:.3f} mS/m")

    print("====================================================\n")
