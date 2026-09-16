# Mixing Rules

This directory contains effective-medium models based on the Bruggeman formulation for estimating the electrical conductivity of cortical and trabecular bone at 100 kHz.

## Cortical Bone

The `bruggeman_cortical.py` script implements a four-phase anisotropic Bruggeman model for cortical bone.

The model considers the cortical bone matrix and vascular channels with different orientations.

To use the model, set the desired bone volume fraction (BV/TV) in the `EXAMPLE` section:

```python
BVTV = 0.90
```

The script returns the effective electrical conductivity in the x, y, and z directions in mS/m.

## Trabecular Bone

The `bruggeman_trabecular.py` script implements a three-phase Bruggeman model for trabecular bone.

The model considers bone matrix, free water, and bone marrow. The free-water fraction is defined as a function of BV/TV.

To use the model, set the desired bone volume fraction (BV/TV) in the `EXAMPLE` section:

```python
BVTV = 0.44
```

The script returns the effective electrical conductivity in mS/m.