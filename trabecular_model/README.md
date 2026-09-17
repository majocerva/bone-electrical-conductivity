# Trabecular Bone Model

This directory contains a representative computational model of trabecular bone and the code used to calculate its effective electrical conductivity.

The model represents a three-dimensional cylindrical trabecular bone sample based on a micro-CT reconstructed geometry. The sample consists of two material domains: bone matrix and bone marrow.

The electrical response is computed using the finite element method (FEM) implemented in FEniCSx.

## Files

- `solver_trabecular.py` – Python code used to compute the effective electrical conductivity of trabecular bone.
- `.msh` file – Representative cylindrical trabecular bone mesh in Gmsh format. The mesh geometry is defined in millimeters (mm).

Detailed information about the mathematical formulation, material domains, boundary conditions, units, calculation of the effective conductivity, and software requirements is provided in the header of `solver_trab.py`.

## Usage

The parameters that can be modified by the user are grouped in the `USER INPUT` section at the end of `solver_trab.py`.

Only the parameters in this section need to be modified to run the model with different trabecular bone samples or material properties.

The user can specify:

- `sample_name` – Name of the Gmsh mesh file (`.msh` extension omitted).
- `sample_length` – Length of the cylindrical sample in mm.
- `sample_radius` – Radius of the cylindrical sample in mm.
- `sigma_matrix` – Electrical conductivity of the bone matrix in mS/m.
- `sigma_marrow` – Electrical conductivity of the bone marrow in mS/m.

For example:

```python
sample_name = "modelo_muestra_3c_231"

sample_length = 2.6
sample_radius = 2.1

sigma_matrix = 20.0
sigma_marrow = 301.0
```

The input Gmsh mesh (`.msh`) must be defined in millimeters (mm). The mesh coordinates are internally converted to meters (m) before solving the FEM problem.

To run the simulation:

```bash
python solver_trab.py
```

The script returns the effective electrical conductivity of the trabecular bone sample in mS/m.

The provided example corresponds to a frequency of **100 kHz**.