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

### CasADi & Yop Dependency

The path planning (reference generation) utilized in the `pathplanner.m` function (Example 2) requires [CasADi](https://web.casadi.org/) and the [Yop toolbox](https://github.com/yoptimalcontrol/yop). 

**Important: Modified Toolbox Version** To implement the time-optimal path planning approach, certain core functions in the Yop toolbox were modified to handle the coordinate transformation presented in the paper above mentioned. 

* **Requirement:** It is imperative to use the specific version of the Yop toolbox provided within this repository. 
* **Warning:** Using a standard or different version of Yop will result in errors or incorrect trajectory generation, as it lacks the necessary modifications required for this formulation to solve correctly. 

### Reference Trajectory

The control of the crane in Example 2 consists of following a reference trajectory obtained using the time-optimal control approach presented in:

Filipe Marques Barbosa, and Johan Löfberg. [*Time-optimal control of cranes subject to container height constraints*](https://doi.org/10.23919/ACC53348.2022.9867816). In 2022 American Control Conference (ACC), pp. 3558-3563. IEEE. doi: 10.23919/ACC53348.2022.9867816.

### Updates in the approximation

In a later version, we updated the exponential-cone approximation to address the asymptotic behavior of the probit function, since $\Phi^{-1}(1-\gamma) \to \infty$ as $\gamma \to 0$. This was done by incorporating a logarithmic term into the approximation. As a result, numerical results may differ (the approximation is more accurate for very small risk levels though).



