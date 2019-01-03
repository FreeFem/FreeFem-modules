---
name: "maxwellWaveguideBC"
category: "electromagnetism, academic"
layout: module
---

# maxwellWaveguideBC

This module illustrates how to implement waveguide boundary conditions in Maxwell equations. So far it is done for 2D planar case in Cartezian and Axial (m = 0) cases. Waveguides modes and waveguides port are hardcoded. 

## Problem
The module uses Maxwell equations written for $\mathbf{E]$ field in $\exp(-i \omega t)$ presentation: $curl curl \mathbf{E} = \frac{\omega^2}{c^2} \epsilon(x,y) \mathbf{E}$ (here: $\omega$ - is the frequency, $c$ - the speed of light, $\epsilon$ - dielectric permittivity; the CGS system of units is used) 
This module illustrates how to set boundary condition so that there is two waveguide ports through one of which an input signal enters. That is it solves a simple waveguide scattering problem

## Variational form
Variational form is 
$$\int_V dV ( curl(\mathbf{E}) curl(\mathbf{V}) - \frac{\omega^2}{c^2} \epsilon(x,y) \mathbf{E} \mathbf{V}) + \int_{dV} dS curl (\mathbf{E}) \mathbf{V} = 0$$
For waveguide ports the last term in the left hand side is written as
$$\int_{dS_{port}} curl \mathbf{E} \mathbf{V} =  \sum_{s=1,N_{modes}} C_s \int_{dS} curl \mathbf{E_s} \mathbf{V}$$
Here C_s is the amplitude of mode $s$ and $\mathbf(E_s)$ - the s-th mode's electric field. The $C_s$ in its turn is written as 
$$C_s = \frac{1}{N_s} \int_{S_{port}} dS ( \mathbf{E} \times \mathbf{H_s} - \mathbf{E_s} \times \mathbf{H})$$, here $\mathbf{H}$ and $\mathbf{H_s}$ are current and mode's magnetic field which are found from Maxwell equations, $N_s$ is the norm. 

Finally, the last term of the LHS side can be written in term of the direct product as
$$
\int_{dS_{port}} curl \mathbf{E} \mathbf{V} =  \sum_{s=1,N_{modes}} \frac{1}{N_s} \int_{S_{port}} dS ( \mathbf{E} \times \mathbf{H_s} - \mathbf{E_s} \times \mathbf{H}) \bigotimes \int_{dS} curl \mathbf{E_s} \mathbf{V},
$$
that allows one to implement it by means of FreeFem++ tools without low-level matrix assembly (for the price of performance, of course).

## Algorithms

[Case description]

_[Optional]_
### 2D

[Case description]

### 3D

[Case description]

### Optional

[e.g., Gmsh script]

_[End optional]_

## Validation

[Validation algorithms and results, e.g. convergence curves]

## References

[Links to references, in open-access if possible]

## Authors

[List of authors, optionnaly with a link to their webpage or email address]
