# Optimizing network robustness via Krylov subspaces

This repository contains the code (mostly matlab) and the data for the network robustness optimization methods presented in 
> Optimizing network robustness via Krylov subspaces, by Stefano Massei and Francesco Tudisco, https://arxiv.org/abs/2303.04971

The repository "Tests" contain 6 scripts that replicates the numerical test reported in the manuscript. In particular:

Discrete optimization
test_unweighted_break ---> greedy strategies for the deletion problem on general and transportation networks (Manuscript's section 5.1)
test_unweighted_make  ---> greedy strategies for the addition problem on general and transportation networks (Manuscript's section 5.2)

Continuous optimization
test_weighted_exp_lbfgs ---> addition, tuning, and rewiring via IP with l-bfgs to optimize trace(exp(A)) (Manuscript's section 6)
test_weighted_exp_hessian ---> addition, tuning, and rewiring via IP with Hessian to optimize trace(exp(A)) (Manuscript's section 6)
test_weighted_sinh_lbfgs ---> addition, tuning, and rewiring via IP with l-bfgs to optimize trace(sinh(A)) (Manuscript's section 6)
test_weighted_sinh_hessian ---> addition, tuning, and rewiring via IP with Hessian to optimize trace(sinh(A)) (Manuscript's section 6)

If you use this repository, please cite our work as:

Massei, Stefano, and Francesco Tudisco. "Optimizing network robustness via Krylov subspaces." arXiv preprint arXiv:2303.04971 (2023).
