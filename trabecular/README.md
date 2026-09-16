# Trabecular Bone Model

This directory contains a representative computational model of trabecular bone used to study its effective electrical conductivity.

The trabecular microstructure is inspired by real micro-computed tomography (μCT) images. The model includes two material domains: bone matrix and bone marrow.

The electrical response is computed using the finite element method (FEM) at a frequency of 100 kHz.

The mesh is provided in Gmsh (`.msh`) format. The corresponding FEM implementation is available in the `scripts/` directory of the main repository.
