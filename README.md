# GP-HMAT-MATLAB
Scalable Gaussian Process Regression with Hierarchical Low Rank Matrices


## Table of contents
* [General info](#general-info)
* [Technical Ingredients](#ingredients)
* [File Contents](#contents)

## General info
This paper will be on arXiv (https://arxiv.org/abs/00000).
## Ingredients
The approach consists of two main ingredients: 1) Hierarchical decomposition of a large matrix which requires aggregation (grouping) of degrees of freedom in different levels of hierarchy. This is achieved via an algebraic multigrid technique. 2) Low rank approximation of large off-diagoanl blocks which is achieved via a randomized SVD approach with interpolative decomsposition. 

The implementation includes hierarchical derivative computations for a standard log likelihood for GP training. Linear solves for n=1e6 and n=1e5 GP nodes on a single CPU takes slightly more than a minute and a few seconds respectively.

The method empirically exhibits O(nlog(n)) scalability for sloving Ax=y, which is in suffiecient agreement with the analytical cost estimate. In some cases, empiricall results are more promising than the analytical cost estimate.  

<img src="matrix_self1.png" width="400">  <img src="tree1.png" width="400" style="vertical-align:top"> 

 
## Contents
HMAT folder contains necessary files for performing linear solve Ax=y. A simple MATLAB script, simple_test.m is provided to demonstrate how the main solver back_solve works. 

GP-HMAT folder contains necessary files for computation of (log) likelihood and its derivative. A simple MATLAB script, simple_test.m is provided to demonstrate how the main solver lkl_eval works. 

Large Dataset Example folder contains all necessary files for likelihood optimization and GP regression. A simple MATLAB script MAIN.m is provided to demonstrate GP likelihood optimization and GP regression associated with a large dataset described in the third numerical example of the paper.  For regression on different datasets, only training and testing datasets should be changed. 

______________________
## Contact
For more information please contact the first author via vkeshava@sci.utah.edu.
