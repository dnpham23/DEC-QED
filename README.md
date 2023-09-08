# DEC-QED
A computational toolbox written in Julia for modeling the dynamics of three-dimensional superconducting devices and materials

## Introduction
DEC-QED combines gauge-invariant, flux-based formulation of electrodynamics with Discrete Exterior Calculus to simulate the nonlinear dynamics of multiscale, Josephson-based superconducting devices in microwave environments. For a thorough discussion of the technical details and key results, please check out our paper [[1]](#1).

## Description of the method
In our model, the equations of motion describing the dynamics between light $(\mathbf{A},V)$ and the collective superconducting condensate $\psi(=\sqrt{\rho}e^{i\theta})$ are written for the gauge-invariant, hybridized field $\mathbf{A}'=\mathbf{A-\frac{\hbar}{q}\nabla\theta}$ 

$${\bf\nabla}\\!\times\\!{\bf\nabla}\\!\times\\!{\bf A'} + \mu_0 \epsilon_0\frac{\partial^2{\bf A'}}{\partial t^2} + \frac{\mu_0 q^2}{m}\rho{\bf A'}  -\frac{\mu_0\epsilon_0 q}{2m}\frac{\partial}{\partial t}\nabla\big|{\bf A'}\big|^2 + \frac{\mu_0\epsilon_0\hbar^2}{2mq}\frac{\partial}{\partial t}\nabla\bigg[\frac{\nabla^2(\sqrt{\rho})}{\sqrt{\rho}}\bigg] =  \mu_0{\bf J}_{src},$$

and the density of superconducting condensate, $\rho$, 

$$\frac{\partial\rho}{\partial t}  = {\bf \nabla} \cdot \Bigg[\frac{q}{m}\rho{\bf A'}  - \frac{{\bf J_{src}}}{q} \Bigg] - \frac{\partial\rho_{src}}{\partial t}.$$

For numerical modeling, we employ a dual-mesh construction (primal + dual meshes), on which the continuous equations above are then transformed into equations for the coarse-grained variables $\rho(v)$ and $\Phi(e) = \int_ed{\bf \ell} \cdot {\bf A}'$, where the field $\Phi(e)$ live on the discrete edges of the discretized primal computational mesh, and $\rho(v)$ lives on the vertices of the primal mesh. This procedure is done using Discrete Exterior Calculus (DEC). A schematic of a cubial dual mesh, along with DEC operators used in our formulation is given in the following Figure:

![DEC schematics](/docs/figs/DualMesh_DEC_schematics_4.svg)

## Development status
DEC-QED as a computational toolbox is still in early stage of development. Examples of applications of the method to specific problems are organized into wrapper files. 
The generalization of the toolbox to perform simulations of arbitrary geometries is currently under development.

## Usage 
First, make sure Git is configured on your workstation. Then clone this repository either via SSH key or HTTPS.
At the moment, the DEC-QED toolbox offers two main types of calculations: (1) semi-classical time-dynamical simulations of superconducting devices in a sourced electromagnetic environment and (2) spectral calculations of closed and open systems.

There are two formulations of DEC implemented: (1) DEC with cubical elements and (2) DEC with tetrahedra elements. 

For a calculation using a tetrahedral mesh, the mesh should be constructed using [Gmsh](https://gmsh.info/). The ".msh" output file from Gmsh is then passed into the mesh parser that extracts useful info and builds the primal and dual mesh structs needed for DEC calculations. 
Checkout the wrappers in [simplicial_decqed](https://github.com/dnpham23/DEC-QED/tree/main/simplicial_decqed) for examples on calculations with simplicial (tetrahedra in 3D) meshes.

For cubical-DEC the construction of the primal and dual meshes needs to be constructed within the program in a case-specific manner. 

If you have questions and/or suggestions regarding the toolbox or its theoretical formulation, please direct them to dnpham@princeton.edu

## References
<a id="1">[1]</a> 
D. N. Pham, W. Fan, M. G. Scheer, and H. E. Tureci, [Phys. Rev. A 107, 053704](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.107.053704) (2023),
arXiv version: [arXiv:2212.12775](https://arxiv.org/abs/2212.12775) (2022).

