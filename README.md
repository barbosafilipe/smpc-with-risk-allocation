# Exponential cone approach to joint chance constraints in stochastic model predictive control

This repository contains code to reproduce the results from the paper:

Filipe Marques Barbosa and Johan Löfberg. [*Exponential Cone Approach to Joint Chance Constraints in Stochastic Model Predictive Control.*](https://doi.org/10.1080/00207179.2025.2492305) International Journal of Control, 2025. doi:10.1080/00207179.2025.2492305.

## Citation

If you use this code or the results in academic work, please cite the paper above.

@article{Barbosa2025,<br>
	title = {Exponential cone approach to joint chance constraints in stochastic model predictive control},<br>
	journal = {International Journal of Control},<br>
	publisher = {Taylor \& Francis},<br>
	author = {Barbosa,  Filipe M. and L{\"{o}}fberg, Johan},<br>
	year = {2025},<br>
	month = Apr,<br>
	pages = {3024–3034},<br>
	volume = {98},<br>
	number = {12}<br>
}

## Information

The main contribution of this paper is a conic convex formulation of chance constraints in which the risk allocation is treated as an optimization variable. This formulation is implemented internally in [YALMIP](https://github.com/yalmip/YALMIP/tree/stochastics_conic).

Since the module for chance constraints has not yet been released in the main branch, the `stochastics_conic` branch must be checked out explicitly.
