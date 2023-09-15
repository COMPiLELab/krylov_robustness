# Optimizing network robustness via Krylov subspaces

This repository contains the code (mostly matlab) and the data for the network robustness optimization methods presented in 
> Optimizing network robustness via Krylov subspaces, by Stefano Massei and Francesco Tudisco, https://arxiv.org/abs/2303.04971

If you use this repository, please cite our work as:

``` 
@article{massei2023optimizing,
  title={Optimizing network robustness via Krylov subspaces},
  author={Massei, Stefano and Tudisco, Francesco},
  journal={arXiv:2303.04971},
  year={2023}
}
```

### HOW TO RUN THE TESTS OF THE MANUSCRIPT

The repository "Tests" contain 8 scripts that replicates the numerical test reported in the manuscript. In particular:

### Discrete optimization
```test_unweighted_break``` &rarr; greedy strategies for the deletion problem on general and transportation networks (Manuscript's section 5.1, Table 2)  
```test_unweighted_break_budget``` &rarr; greedy strategies for the deletion problem on transportation networks with varying budget (Manuscript's section 5.1, Figures 1 and 2)  
```test_unweighted_make```  &rarr; greedy strategies for the addition problem on general and transportation networks (Manuscript's section 5.2, Table 3)  
```test_unweighted_make_budget``` &rarr; greedy strategies for the addition problem on transportation networks with varying budget (Manuscript's section 5.2, Figures 3 and 4)  

### Continuous optimization  
```test_weighted_exp_lbfgs``` &rarr; addition, tuning, and rewiring via IP with l-bfgs to optimize trace(exp(A)) (Manuscript's section 6, Table 5)   
```test_weighted_exp_hessian``` &rarr; addition, tuning, and rewiring via IP with Hessian to optimize trace(exp(A)) (Manuscript's section 6, Table 5)  
```test_weighted_sinh_lbfgs``` &rarr; addition, tuning, and rewiring via IP with l-bfgs to optimize trace(sinh(A)) (Manuscript's section 6, Table 6)  
```test_weighted_sinh_hessian``` &rarr; addition, tuning, and rewiring via IP with Hessian to optimize trace(sinh(A)) (Manuscript's section 6, Table 6)  

### EXTERNAL CODE

The repository "MIOBI Codes" contains the implementation of the greedy algorithms proposed in 

```
Chan, Hau, Leman Akoglu, and Hanghang Tong. "Make it or break it: Manipulating robustness in large networks." Proceedings of the 2014 SIAM International Conference on Data Mining. Society for Industrial and Applied Mathematics, 2014.
```

that has been kindly forwarded to us by the authors via private communication.

The files ```expmv.m```, ```normAm.m```, ```select_taylor_degree.m```, and ```theta_taylor.mat``` contained in the repository "functions" implement the ```expmv``` algorithm proposed in 
```
Al-Mohy, Awad H., and Nicholas J. Higham. "Computing the action of the matrix exponential, with an application to exponential integrators." SIAM journal on scientific computing 33.2 (2011): 488-511.
```
that is freely available at [https://github.com/higham/expmv](https://github.com/higham/expmv).
