# Cortical Bone Model

This directory contains a representative computational model of cortical bone used to study its effective electrical conductivity.

The cortical bone geometry is parametrically generated and represents the main components of the bone vascular system, including Haversian and Volkmann canals.

The electrical response is computed using the finite element method (FEM) at a frequency of 100 kHz in the axial and transverse directions.

The mesh is provided in Gmsh (`.msh`) format. The corresponding FEM implementation is available in the `scripts/` directory of the main repository.
