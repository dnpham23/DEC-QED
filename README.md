# DEC-QED
A computational toolbox for modeling the dynamics of three-dimensional superconducting devices and materials

## Introduction
DEC-QED combines gauge-invariant, flux-based formulation of electrodynamics with Discrete Exterior Calculus to simulate the nonlinear dynamics of multiscale, Josephson-based superconducting devices in microwave environments. For a thorough discussion of the technical details and key results, please check out our paper: DN Pham, W Fan, MG Scheer, and HE Tureci, [arXiv:2212.12775](https://arxiv.org/abs/2212.12775) (2022).

## Description
In our model, the dynamical equations describing light-matter interaction are written for the gauge-invariant, hybridized field $\mathbf{A}'=\mathbf{A-\frac{\hbar}{q}\nabla\theta}$ and the density of superconducting condensate, $\rho$. 

$${\bf\nabla}\\!\times\\!{\bf\nabla}\\!\times\\!{\bf A'} + \mu_0 \epsilon_0\frac{\partial^2{\bf A'}}{\partial t^2} + \frac{\mu_0 q^2}{m}\rho{\bf A'}  -\frac{\mu_0\epsilon_0 q}{2m}\frac{\partial}{\partial t}\nabla\big|{\bf A'}\big|^2 + \frac{\mu_0\epsilon_0\hbar^2}{2mq}\frac{\partial}{\partial t}\nabla\bigg[\frac{\nabla^2(\sqrt{\rho})}{\sqrt{\rho}}\bigg] =  \mu_0{\bf J}_{src},$$

$$\frac{\partial\rho}{\partial t}  = {\bf \nabla} \cdot \Bigg[\frac{q}{m}\rho{\bf A'} - \frac{{\bf J}_{src}}{q}\Bigg] - \frac{\partial\rho_{src}}{\partial t}.$$

For numerical modeling, we employ a dual-mesh construction (primal + dual meshes), on which the continuous equations above are then transformed into equations for the coarse-grained variables $\rho(v)$ and $\Phi(e) = \int_ed{\bf \ell} \cdot {\bf A}'$, where the field $\Phi(e)$ live on the discrete edges of the discretized primal computational mesh, and $\rho(v)$ lives on the vertices of the primal mesh. This procedure is done using Discrete Exterior Calculus (DEC) [[1]](#1). A schematic of a cubial dual mesh, along with DEC operators used in our formulation is given in the following Figure:

## References
<a id="1">[1]</a> 
Hirani, A. N.. (2003). 
Go to statement considered harmful. 
PhD Thesis.
