# Bone Electrical Conductivity

This repository contains supplementary data, representative computational models, and scripts developed to study the electrical properties of bone tissue.

The repository is organized into trabecular, cortical, mixing-rule, and macroscale models.

## Trabecular Bone

Computational models of trabecular bone with microstructures inspired by real micro-computed tomography (μCT) images. The electrical response is obtained using the finite element method (FEM).

A representative trabecular model and the corresponding computational scripts are provided.

## Cortical Bone

Computational models of cortical bone with parametric geometries representing the main components of the bone vascular system. The electrical response is obtained using the finite element method (FEM) in the axial and transverse directions.

A representative cortical model and the corresponding computational scripts are provided.

## Mixing Rules

Effective medium approaches based on the Bruggeman formulation are included to estimate the effective electrical conductivity of bone tissue and complement the FEM results.

## Macroscale

A macroscale computational bone model integrates the electrical properties of cortical bone, trabecular bone, and bone marrow. The electrical response is obtained using the finite element method (FEM) considering different tissue conditions.

## Supplementary Data

The numerical results used in the manuscript are provided in:

`Supplementary_Data.xlsx`

The workbook includes FEM results and literature data for trabecular and cortical bone.

## Repository Structure

- `trabecular/` – Representative trabecular model and FEM scripts.
- `cortical/` – Representative cortical model and FEM scripts.
- `mixing_rules/` – Effective medium models and associated scripts.
- `macroscale/` – Macroscale model and associated scripts.
- `Supplementary_Data.xlsx` – Numerical data associated with the manuscript.

## Frequency

The electrical results presented in this repository correspond to a frequency of **100 kHz**.
