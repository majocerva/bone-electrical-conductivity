# Macroscale Tibia Model

This directory contains a representative computational model of the tibia and the code used to calculate its electrical impedance.

The model represents a three-dimensional tibia geometry composed of three material domains: cortical bone, trabecular bone, and bone marrow.

The electrical response is computed using the finite element method (FEM) implemented in FEniCSx. Cortical bone is represented using an anisotropic conductivity tensor, whereas trabecular bone and bone marrow are considered isotropic.

## Files

- `solver_tibia.py` – Python code used to compute the electrical impedance of the tibia model.
- `.msh` file – Representative three-dimensional tibia mesh in Gmsh format. The mesh geometry is defined in meters (m).

Detailed information about the mathematical formulation, material domains, boundary conditions, units, calculation of the electrical impedance, and software requirements is provided in the header of `solver_tibia.py`.

## Usage

The parameters that can be modified by the user are grouped in the `USER INPUT` section at the end of `solver_tibia.py`.

Only the parameters in this section need to be modified to run the model with different tibia geometries, electrode configurations, or material properties.

The user can specify:

- `sample_name` – Name of the Gmsh mesh file (`.msh` extension omitted).
- `configuration` – Electrode configuration: `"H"` for horizontal or `"V"` for vertical.
- `sigma_cortical_x` – Cortical bone conductivity in the transverse direction in mS/m.
- `sigma_cortical_z` – Cortical bone conductivity in the longitudinal direction in mS/m.
- `sigma_trabecular` – Trabecular bone conductivity in mS/m.
- `sigma_marrow` – Bone marrow conductivity in mS/m.

For example:

```python
sample_name = "Tibia_SANA_5mm_3c"

configuration = "H"

sigma_cortical_x = 6.363
sigma_cortical_z = 11.062

sigma_trabecular = 120.395
sigma_marrow = 300.0
```

The input Gmsh mesh (`.msh`) must be defined in meters (m).

To run the simulation:

```bash
python solver_macro.py
```

The script returns the electrical impedance of the tibia model in kΩ and the active and passive electrode areas in mm².

The conductivity values used in the provided example correspond to a frequency of **100 kHz**.