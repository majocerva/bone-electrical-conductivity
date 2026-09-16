# Trabecular Bone

This directory contains a representative computational model of trabecular bone used to study its effective electrical conductivity.

The trabecular microstructure is inspired by real micro-computed tomography (μCT) images. The electrical response is obtained using the finite element method (FEM) at a frequency of 100 kHz.

## Model

The `models/` directory contains a representative trabecular bone geometry used in the FEM simulations.

The model consists of two material domains:

- Bone matrix
- Bone marrow

The corresponding electrical properties and physical tags are defined in the simulation scripts available in the `scripts/` directory of the main repository.
