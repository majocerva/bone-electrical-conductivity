# Bone Electrical Conductivity

This repository contains numerical models, data, and scripts developed to study the electrical properties of bone tissue using a multiscale modeling approach.

The work is organized into three modeling levels: trabecular, cortical, and macroscale.

## Trabecular Bone

Computational models of trabecular bone with microstructures inspired by real micro-computed tomography (μCT) images. The electrical response is obtained using the finite element method (FEM).

The numerical results are complemented by an effective medium model based on the Bruggeman formulation to analyze the relationship between bone microstructure and effective electrical conductivity.

## Cortical Bone

Computational models of cortical bone with parametric geometries representing the main components of the bone vascular system. The electrical response is obtained using the finite element method (FEM), allowing the relationship between porosity and electrical conductivity in the axial and transverse directions to be investigated.

The numerical results are complemented by an effective medium model based on the Bruggeman formulation.

## Macroscale

A macroscale computational bone model integrating the electrical properties of cortical bone, trabecular bone, and bone marrow. The electrical response is obtained using the finite element method (FEM), considering different tissue conditions.

## Repository Structure

- `trabecular/` – Models, data, and scripts related to trabecular bone.
- `cortical/` – Models, data, and scripts related to cortical bone.
- `macroscale/` – Models, data, and scripts related to the macroscale model.

## Frequency

The electrical results presented in this repository correspond to a frequency of **100 kHz**.
