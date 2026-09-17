# Cortical Bone FEM Model

This directory contains a representative computational model of cortical bone and the code used to calculate its effective electrical conductivity.

The model represents a three-dimensional cubic cortical bone microstructure with a parametrically generated geometry. The model consists of cortical bone matrix and vascular channels, including Haversian and Volkmann canals.

The electrical response is computed using the finite element method (FEM) implemented in FEniCSx.

## Files

- `solver_cortical.py` – Python code used to compute the effective electrical conductivity of cortical bone.
- `.msh` file – Representative cubic cortical bone mesh in Gmsh format. The mesh geometry is defined in meters (m).

Detailed information about the mathematical formulation, material domains, boundary conditions, units, calculation of the effective conductivity, and software requirements is provided in the header of `solver_cortical.py`.

## Usage

The parameters that can be modified by the user are grouped in the `USER INPUT` section at the end of `solver_cortical.py`.

Only the parameters in this section need to be modified to run the model with different cortical bone samples, material properties, or electric field directions.

The user can specify:

- `sample_name` – Name of the Gmsh mesh file (`.msh` extension omitted).
- `direction` – Applied electric field direction (`x`, `y`, or `z`).
- `sample_length` – Length of the cubic sample in m.
- `sigma_matrix` – Electrical conductivity of the cortical bone matrix in mS/m.
- `sigma_marrow` – Electrical conductivity of the vascular channels in mS/m.

For example:

```python
sample_name = "Poro6"

direction = "y"

sample_length = 1.0e-3

sigma_matrix = 3.84
sigma_marrow = 267.0
```

The input Gmsh mesh (`.msh`) must be defined in meters (m).

To run the simulation:

```bash
python solver_cortical.py
```

The script returns the effective electrical conductivity of the cortical bone model in mS/m.

The provided example corresponds to a frequency of **100 kHz**.