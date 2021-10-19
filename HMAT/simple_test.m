%%%%%For full documentation of the approach, please see the Arxiv
%%%%%paper: GP-HMAT: SCALABLE, O(n log(n)) GAUSSIAN PROCESS
%%%%%REGRESSION WITH HIERARCHICAL LOW RANK MATRICES 

clc;clear all;close all;

warning('off')

load data;
%%nodes are generated with rand
%%y is also generated with rand and normalized to 1.

n=5e3;nodes_test=nodes(:,1:n);y_test=y(1:n);y_test=y_test/norm(y_test);

%%nkernel determines exponetial kernel (1) or squared exponential kernel (2)
%l is the hyperparameter; either scalar value or column vector with size
%equal to data dimension
%nmode detrmines computation of derivative; in HMAT, i.e. back_solve function, it is set as 0.
%if you need derivative of the solution with respect to hyperparameters see
%lkl_eval function in GP-HMAT
nkernel=2;l=[1];mode=0;%you may change l to l=[1;1] to test for ARD kernel 

kernel_params={nkernel l mode}; %parameter of the kernel
cutoff_size=1005;%cut off size for matrix slicing
k=30;%rank parameter
delta1=1e-3; %regularization for diagonal matrices
delta2=0; %regularization for SMW correction matrix; always zero recommended

params={kernel_params,cutoff_size,k,delta1,delta2};


%%%%%%HMAT Computations%%%%%%%%%%%%%%%%%%%%%
tic;
[sol_HMAT]=back_solve(nodes_test,{y_test},params);
t_HMAT=toc;
fprintf('The HMAT solver time is %s.\n',t_HMAT);
%%%%%%%%%%%%%%%%%%%%%%%MATLAB backslash \ Computations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%changing kernel params to kernel_params={nkernel l 1} provides matrix
%derivative in case needed. This only works in kernel_computation.
%if you need derivative of the solution with respect to hyperparameters see
%lkl_eval function in GP-HMAT

tic;
[K_MATLAB]=kernel_computation(nodes_test,nodes_test,kernel_params);
sol_MATLAB=(K_MATLAB+delta1*eye(n))\y_test;
t_MATLAB=toc;
fprintf('The MATLAB solver time is %s.\n',t_MATLAB)

%%%%%%%Wall-clock time
%t_HMAT=0.2809
%t_MATLAB=1.0695

%%%%Normalized Error
norm_err=norm(sol_HMAT-sol_MATLAB)/norm(sol_MATLAB);
fprintf('The normalized error is %s.\n',norm_err)
%%%in the order of 1e-8



