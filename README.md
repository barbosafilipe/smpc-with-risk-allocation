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

## Information for the user

The main contribution of this paper is the formulation of chance constraints as (exponential) conic convex functions, which allows risk allocation to be treated as an optimization variable.

### YALMIP dependency

Once chance constraints are added using `probability()`, the formulation is handled internally by [YALMIP](https://github.com/yalmip/YALMIP/tree/stochastics_conic). Therefore, YALMIP must be installed and added to your path.

**Note:** The module for chance-constrained optimization has not yet been officially released in the main branch. The modeling is implemented in the `stochastics_conic` branch, which must be used.

### MOSEK dependency

The solver used is [MOSEK](https://www.mosek.com/), which provides free academic [licenses](https://www.mosek.com/products/academic-licenses/).

### Updates in the approximation

In a later version, we updated the exponential-cone approximation to address the asymptotic behavior of the probit function, since $\Phi^{-1}(1-\gamma) \to \infty$ as $\gamma \to 0$. This was done by incorporating a logarithmic term into the approximation. As a result, numerical results may differ (the approximation is more accurate for very small risk levels though).



