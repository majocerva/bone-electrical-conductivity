# Cortical Bone Model

This directory contains a representative computational model of cortical bone, the code used to generate the parametric geometry and finite element mesh, and the code used to calculate its effective electrical conductivity.

The model represents a three-dimensional cubic cortical bone microstructure with a parametrically generated geometry. The model consists of cortical bone matrix and vascular channels, including Haversian and Volkmann canals.

The geometry and mesh are generated using Gmsh, while the electrical response is computed using the finite element method (FEM) implemented in FEniCSx.

## Files

* `cortical.py` – Original code containing the functions used to generate the parametric cortical bone geometry, calculate porosity, and generate the Gmsh mesh.
* `generate_mesh.py` – User input script used to define the geometric parameters and generate a cortical bone geometry and mesh.
* `solver_cortical.py` – Python code used to compute the effective electrical conductivity of cortical bone.
* `.msh` file – Representative cubic cortical bone mesh in Gmsh format. The mesh geometry is defined in meters (m).

Detailed information about the mathematical formulation, material domains, boundary conditions, units, calculation of the effective conductivity, and software requirements is provided in the header of `solver_cortical.py`.

## Mesh Generation

A new cortical bone geometry and mesh can be generated using `generate_mesh.py`.

The geometric parameters that can be modified by the user are grouped in the `USER INPUT` section of `generate_mesh.py`. These include:

* Sample size and osteon density.
* Minimum and maximum osteon diameters.
* Minimum and maximum Haversian canal diameters and inclination angles.
* Minimum and maximum Volkmann canal diameters and inclination angles.
* Random seed used for the parametric geometry generation.
* Sample name.

The script calls the geometry and mesh generation functions defined in `cortical.py`. Therefore, `generate_mesh.py` and `cortical_geometry.py` must be located in the same directory.

To generate a new geometry and mesh:

```bash
python generate_mesh.py
```

The script generates:

* a `.geo` file containing the parametric Gmsh geometry,
* a `.msh` file containing the finite element mesh,
* the porosity of the generated cortical bone sample.

The generated `.msh` file is scaled to meters (m) and can be used directly as input for `solver_cortical.py`.

A representative `.msh` file is also provided in this directory, so mesh generation is not required to reproduce the example FEM simulation.

## FEM Simulation

The parameters that can be modified by the user are grouped in the `USER INPUT` section at the end of `solver_cortical.py`.

Only the parameters in this section need to be modified to run the model with different cortical bone samples, material properties, or electric field directions.

The user can specify:

* `sample_name` – Name of the Gmsh mesh file (`.msh` extension omitted).
* `direction` – Applied electric field direction (`x`, `y`, or `z`).
* `sample_length` – Length of the cubic sample in m.
* `sigma_matrix` – Electrical conductivity of the cortical bone matrix in mS/m.
* `sigma_marrow` – Electrical conductivity of the vascular channels in mS/m.

For example:

```python
sample_name = "Poro6"

direction = "y"

sample_length = 1.0e-3

sigma_matrix = 3.84
sigma_marrow = 267.0
```

The input Gmsh mesh (`.msh`) must be defined in meters (m).

To run the FEM simulation:

```bash
python solver_cortical.py
```

The script returns the effective electrical conductivity of the cortical bone model in mS/m.

The provided example corresponds to a frequency of **100 kHz**.
