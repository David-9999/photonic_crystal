# Photonic_crystal

This program can compute intensity of electirc field and transmition coeficient of 1D photonic crystal with defect.

## Theory


Theory is based on $vloz link na knihu$.
In simplicity photonic crystal is material that disrupts translation invariance of free space and the index of refraction is periodically repeated. The simpliest photonic crystal is repeatiog of two layers of materials. We call this 1D photonic crystal.

Defect is everything that disrupts periodically repeated of index of refraction. For example photonic crystal with defect can be ilustrated by N layers of 2 materials (A and B), then D layer and another N layers of A and B.

We gonna need some physical expressions from electromagnetics and optics.
The transition of electric field can be described with matrix of transition.

$$ {\left\lbrack \matrix{E_2^+ \cr E_2^-} \right\rbrack}=M_{\text{T}}{\left\lbrack \matrix{E_1^+ \cr E_1^-} \right\rbrack} $$

The matrix of transition can be created by combination of two types of matries. 
First is matrix that describe transmition through layer,

$$ M_{\text{layer}}={\left\lbrack \matrix{e^{ik_z\Delta} & 0 \cr 0 & e^{-ik_z\Delta}} \right\rbrack} $$

and maxtrix of transmition throght interfrace of layers,

$$ M_{\text{TE}}={\left\lbrack \matrix{1+χ_E & 1-χ_E \cr 1-χ_E & 1+χ_E} \right\rbrack} $$

where
$$\chi_E = \frac{\mu_2 k_{1z}}{\mu_1 k_{2z}}$$
$$k_{1,2z}=\frac{\omega}{c}n_{1,2}\cos\left(\theta_{1,2}\right)$$

$\omega$ is angle of incidence, $\mu$ is permeability and $n=\sqrt{\epsilon \mu}$ is refraction index 
First of all you need to input permitivity, permeability, lenghts of layers, angle of incidence, size of defect and number of repeating of A+B layer.
Calculations are based on


