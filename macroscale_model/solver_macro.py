#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 10:03:22 2026

@author: mjcervantes
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
===============================================================================
FEniCSx – Electrical Impedance of a Multiscale Tibia Model
===============================================================================

This script computes the electrical impedance of a three-dimensional in silico
tibia model composed of cortical bone, trabecular bone, and marrow.

The model is represented as a heterogeneous three-phase domain composed of:

    - Cortical bone
    - Trabecular bone
    - Marrow

The electrical potential distribution is obtained by solving the steady-state
conduction equation

        div(sigma * grad(V)) = 0

using the finite element method (FEM) implemented in FEniCSx.

The cortical bone conductivity is represented by an anisotropic diagonal
conductivity tensor, whereas trabecular bone and marrow are considered
isotropic.

Dirichlet boundary conditions are imposed on two surface electrodes identified
by physical facet tags in the Gmsh mesh.

Two electrode configurations can be selected:

    - Horizontal configuration (H): electrode tags 10 and 20
    - Vertical configuration   (V): electrode tags 30 and 40

The electrical impedance is calculated from the total electrical power
dissipated in the domain according to

        Z = DeltaV^2 / P

where DeltaV is the imposed potential difference and P is the total electrical
power dissipated in the domain.

The tibia geometry is defined in meters, and all FEM calculations are performed
using SI units.

Units
-----
Input geometry                  : m
Internal FEM geometry           : m
Input conductivity              : mS/m
Internal FEM conductivity       : S/m
Impedance output                : kOhm
Electrode area output           : mm²

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

from mpi4py import MPI
from petsc4py import PETSc

from dolfinx import fem
from dolfinx.io import gmshio
from dolfinx.fem.petsc import LinearProblem


# =============================================================================
# CONDUCTIVITY TENSOR
# =============================================================================

def conductivity_tensor(sigma_diag):
    """
    Create a diagonal conductivity tensor.

    Parameters
    ----------
    sigma_diag : array-like
        Conductivity components [sigma_x, sigma_y, sigma_z] in mS/m.

    Returns
    -------
    numpy.ndarray
        3 x 3 conductivity tensor in S/m.
    """

    sigma_diag = np.asarray(sigma_diag, dtype=float) * 1e-3

    tensor = np.zeros((3, 3), dtype=float)

    tensor[0, 0] = sigma_diag[0]
    tensor[1, 1] = sigma_diag[1]
    tensor[2, 2] = sigma_diag[2]

    return tensor


# =============================================================================
# FEM SOLVER
# =============================================================================

def solve_tibia_model(
        name,
        configuration,
        sigma_cortical_x,
        sigma_cortical_z,
        sigma_trabecular,
        sigma_marrow):
    """
    Compute the electrical impedance of the multiscale tibia model.

    Parameters
    ----------
    name : str
        Name of the Gmsh mesh file without the ".msh" extension.

    configuration : str
        Electrode configuration: "H" for horizontal or "V" for vertical.

    sigma_cortical_x : float
        Cortical bone conductivity in the transverse direction in mS/m.

    sigma_cortical_z : float
        Cortical bone conductivity in the longitudinal direction in mS/m.

    sigma_trabecular : float
        Trabecular bone conductivity in mS/m.

    sigma_marrow : float
        Marrow conductivity in mS/m.

    Returns
    -------
    impedance : float
        Electrical impedance in kOhm.

    electrode_areas : tuple
        Active and passive electrode areas in mm².
    """


    # -------------------------------------------------------------------------
    # Geometry
    # -------------------------------------------------------------------------

    msh_file = os.path.join(f"{name}.msh")


    # -------------------------------------------------------------------------
    # Applied potential
    # -------------------------------------------------------------------------

    U0 = 40.0
    Ug = 0.0

    delta_V = U0 - Ug


    # -------------------------------------------------------------------------
    # Mesh
    # -------------------------------------------------------------------------

    msh, cell_tags, facet_tags = gmshio.read_from_msh(
        msh_file,
        MPI.COMM_WORLD,
        0,
        gdim=3
    )

    tdim = msh.topology.dim
    fdim = tdim - 1

    msh.topology.create_connectivity(fdim, tdim)

    V = fem.functionspace(msh, ("P", 1))


    # -------------------------------------------------------------------------
    # Selected electrode configuration
    # -------------------------------------------------------------------------

    configuration = configuration.upper()

    electrode_tags = {
        "H": (10, 20),
        "V": (30, 40)
    }

    if configuration not in electrode_tags:
        raise ValueError(
            "configuration must be 'H' (horizontal) or 'V' (vertical)"
        )

    active_tag, passive_tag = electrode_tags[configuration]


    # -------------------------------------------------------------------------
    # Boundary conditions
    # -------------------------------------------------------------------------

    active_facets = facet_tags.find(active_tag)
    passive_facets = facet_tags.find(passive_tag)

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
    # Electrode areas
    # -------------------------------------------------------------------------

    ds = ufl.Measure(
        "ds",
        domain=msh,
        subdomain_data=facet_tags
    )

    one = fem.Constant(
        msh,
        PETSc.ScalarType(1.0)
    )

    area_active_local = fem.assemble_scalar(
        fem.form(
            one * ds(active_tag)
        )
    )

    area_passive_local = fem.assemble_scalar(
        fem.form(
            one * ds(passive_tag)
        )
    )

    area_active = msh.comm.allreduce(
        area_active_local,
        op=MPI.SUM
    )

    area_passive = msh.comm.allreduce(
        area_passive_local,
        op=MPI.SUM
    )


    # -------------------------------------------------------------------------
    # Material conductivities
    # -------------------------------------------------------------------------

    # Cortical bone: anisotropic
    cortical = conductivity_tensor(
        [
            sigma_cortical_x,
            sigma_cortical_x,
            sigma_cortical_z
        ]
    )

    # Trabecular bone: isotropic
    trabecular = conductivity_tensor(
        [
            sigma_trabecular,
            sigma_trabecular,
            sigma_trabecular
        ]
    )

    # Marrow: isotropic
    marrow = conductivity_tensor(
        [
            sigma_marrow,
            sigma_marrow,
            sigma_marrow
        ]
    )


    # -------------------------------------------------------------------------
    # Spatial conductivity distribution
    # -------------------------------------------------------------------------

    Q = fem.functionspace(
        msh,
        ("DG", 0, (3, 3))
    )

    sigma = fem.Function(
        Q,
        dtype=np.float64
    )


    # Physical tags defined in the Gmsh mesh
    CORTICAL_TAG = 100
    TRABECULAR_TAG = 200
    MARROW_TAG = 300


    cortical_cells = cell_tags.find(
        CORTICAL_TAG
    )

    trabecular_cells = cell_tags.find(
        TRABECULAR_TAG
    )

    marrow_cells = cell_tags.find(
        MARROW_TAG
    )


    def assign_tensor(function, cells, tensor):
        """
        Assign a conductivity tensor to a set of mesh cells.
        """

        values = np.asarray(
            tensor,
            dtype=PETSc.ScalarType
        ).reshape(9)

        dofmap = function.function_space.dofmap

        for cell in cells:

            dofs = dofmap.cell_dofs(cell)

            function.x.array[
                dofs[0] * 9:(dofs[0] + 1) * 9
            ] = values


    assign_tensor(
        sigma,
        cortical_cells,
        cortical
    )

    assign_tensor(
        sigma,
        trabecular_cells,
        trabecular
    )

    assign_tensor(
        sigma,
        marrow_cells,
        marrow
    )


    # -------------------------------------------------------------------------
    # Variational problem
    # -------------------------------------------------------------------------

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    f = fem.Constant(
        msh,
        PETSc.ScalarType(0)
    )


    a = (
        ufl.inner(
            ufl.dot(
                sigma,
                ufl.grad(u)
            ),
            ufl.grad(v)
        )
        * ufl.dx
    )


    L = (
        ufl.inner(
            f,
            v
        )
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
            "ksp_type": "gmres",
            "pc_type": "ilu",
            "ksp_rtol": 1e-10,
            "ksp_atol": 1e-15,
        },
    )


    potential = problem.solve()


    # -------------------------------------------------------------------------
    # Electrical power
    # -------------------------------------------------------------------------

    power_local = fem.assemble_scalar(
        fem.form(
            ufl.inner(
                sigma * ufl.grad(potential),
                ufl.grad(potential)
            )
            * ufl.dx
        )
    )

    power = msh.comm.allreduce(
        power_local,
        op=MPI.SUM
    )


    # -------------------------------------------------------------------------
    # Electrical impedance
    # -------------------------------------------------------------------------

    impedance = delta_V**2 / power

    # Convert Ohm to kOhm
    impedance_kohm = impedance / 1000.0

    # Convert electrode areas from m² to mm²
    area_active_mm2 = area_active * 1e6
    area_passive_mm2 = area_passive * 1e6

    electrode_areas = (
        area_active_mm2,
        area_passive_mm2
    )


    return impedance_kohm, electrode_areas


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
    sample_name = "Tibia_SANA_5mm_3c"


    # -------------------------------------------------------------------------
    # Electrode configuration
    # -------------------------------------------------------------------------

    # Electrode configuration: "H" (horizontal) or "V" (vertical)
    configuration = "H"


    # -------------------------------------------------------------------------
    # Material conductivities [mS/m]
    # -------------------------------------------------------------------------

    # Cortical bone (anisotropic)
    sigma_cortical_x = 6.363
    sigma_cortical_z = 11.062

    # Trabecular bone (isotropic)
    sigma_trabecular = 120.395 

    # Marrow (isotropic)
    sigma_marrow = 300.0



    # =========================================================================
    # RUN FEM MODEL
    # =========================================================================

    impedance, electrode_areas = solve_tibia_model(
        name=sample_name,
        configuration=configuration,
        sigma_cortical_x=sigma_cortical_x,
        sigma_cortical_z=sigma_cortical_z,
        sigma_trabecular=sigma_trabecular,
        sigma_marrow=sigma_marrow
    )


    # =========================================================================
    # RESULTS
    # =========================================================================

    print("\n====================================================")
    print("MULTISCALE TIBIA FEM RESULTS")
    print("====================================================")

    print(f"Sample                    : {sample_name}")
    print(f"Electrode configuration   : {configuration.upper()}")

    print("----------------------------------------------------")

    print(f"Cortical conductivity X   : {sigma_cortical_x:.3f} mS/m")
    print(f"Cortical conductivity Z   : {sigma_cortical_z:.3f} mS/m")
    print(f"Trabecular conductivity   : {sigma_trabecular:.3f} mS/m")
    print(f"Marrow conductivity       : {sigma_marrow:.3f} mS/m")

    print("----------------------------------------------------")

    print(f"Active electrode area     : {electrode_areas[0]:.3f} mm²")
    print(f"Passive electrode area    : {electrode_areas[1]:.3f} mm²")

    print("----------------------------------------------------")

    print(f"Electrical impedance      : {impedance:.3f} kOhm")

    print("====================================================\n")