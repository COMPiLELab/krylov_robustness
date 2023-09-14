# Optimizing network robustness via Krylov subspaces

This repository contains the code (mostly matlab) and the data for the network robustness optimization methods presented in 
> Optimizing network robustness via Krylov subspaces, by Stefano Massei and Francesco Tudisco, https://arxiv.org/abs/2303.04971

If you use this repository, please cite our work as:

Massei, Stefano, and Francesco Tudisco. "Optimizing network robustness via Krylov subspaces." arXiv preprint arXiv:2303.04971 (2023).

HOW TO RUN THE TESTS OF THE MANUSCRIPT

The repository "Tests" contain 6 scripts that replicates the numerical test reported in the manuscript. In particular:

Discrete optimization
test_unweighted_break ---> greedy strategies for the deletion problem on general and transportation networks (Manuscript's section 5.1)
test_unweighted_make  ---> greedy strategies for the addition problem on general and transportation networks (Manuscript's section 5.2)

Continuous optimization
test_weighted_exp_lbfgs ---> addition, tuning, and rewiring via IP with l-bfgs to optimize trace(exp(A)) (Manuscript's section 6)
test_weighted_exp_hessian ---> addition, tuning, and rewiring via IP with Hessian to optimize trace(exp(A)) (Manuscript's section 6)
test_weighted_sinh_lbfgs ---> addition, tuning, and rewiring via IP with l-bfgs to optimize trace(sinh(A)) (Manuscript's section 6)
test_weighted_sinh_hessian ---> addition, tuning, and rewiring via IP with Hessian to optimize trace(sinh(A)) (Manuscript's section 6)

EXTERNAL CODE

The repository "MIOBI Codes" contains the implementation of the greedy algorithms proposed in 

Chan, Hau, Leman Akoglu, and Hanghang Tong. "Make it or break it: Manipulating robustness in large networks." Proceedings of the 2014 SIAM International Conference on Data Mining. Society for Industrial and Applied Mathematics, 2014.

that has kindly forwarded to us by the authors in a private communication.

The files expmv.m, normAm.m, select_taylor_degree.m, and theta_taylor.mat contained in the repository "functions" implement the expmv algorithm proposed in 

Al-Mohy, Awad H., and Nicholas J. Higham. "Computing the action of the matrix exponential, with an application to exponential integrators." SIAM journal on scientific computing 33.2 (2011): 488-511.

that is freely available at https://github.com/higham/expmv.
