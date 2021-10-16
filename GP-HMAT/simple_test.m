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
%nmode detrmines computation of derivative; in GP-HMAT, i.e. lkl_eval function, it is set as 1.
%The derivatives of energy and determinant are already computed in
%energy_d and xdet_d. If you need derivative of the solution with respect to hyperparameters see
%Remark 2 below. 

nkernel=2;l=[1];mode=1; %%you may change l to l=[1;1] to test for ARD kernel 

kernel_params={nkernel l mode}; %parameter of the kernel
cutoff_size=1005;%cut off size for matrix slicing
k=25;%rank parameter
delta1=1e-3; %regularization for diagonal matrices
delta2=0; %regularization for SMW correction matrix; always zero recommended

params={kernel_params,cutoff_size,k,delta1,delta2};

%%%%%%GP-HMAT Computations%%%%%%%%%%%%%%%%%%%%%
tic;
[xdet,xdet_d,energy,energy_d,sol_out,A_d_sol]=lkl_eval(nodes_test,{y_test},params);
t_GPHMAT=toc;
fprintf('The GP-HMAT solver time is %s.\n',t_GPHMAT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%MATLAB Energy and Determinant Computations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
[K_MATLAB]=kernel_computation(nodes_test,nodes_test,kernel_params);
sol_MATLAB=(K_MATLAB{1}+delta1*eye(n))\y_test;
energy_MATLAB=sol_MATLAB'*y_test;
det_MATLAB=sum(log(eig(K_MATLAB{1}+delta1*eye(n))));
%%Derivative business
det_MATLAB_d=zeros(1,length(l));
energy_MATLAB_d=zeros(1,length(l));
A_d_sol_MATLAB=zeros(n,length(l));
for i=1:length(l)
    det_MATLAB_d(1,i)=trace([K_MATLAB{1}+delta1*eye(n)]\K_MATLAB{i+1});
    energy_MATLAB_d(1,i)=-sol_MATLAB'*K_MATLAB{i+1}*sol_MATLAB;
    A_d_sol_MATLAB(:,i)=K_MATLAB{i+1}*sol_MATLAB;
end
t_MATLAB=toc;
fprintf('The MATLAB solver time is %s.\n',t_MATLAB)

%%%%%%%Wall-clock time
%t_GPHMAT=9.15e-1
%t_MATLAB=7.13e+0

%%%%Normalized Error

norm_err_ene=norm(energy-energy_MATLAB)/norm(energy_MATLAB);
norm_err_det=norm(real(xdet)-det_MATLAB)/norm(det_MATLAB);
norm_err_ene_d=norm(energy_d{1}-energy_MATLAB_d(1))/norm(energy_MATLAB_d(1));
norm_err_det_d=norm(xdet_d{1}-det_MATLAB_d(1))/norm(det_MATLAB_d(1));

fprintf('The normalized error in energy is %s.\n', norm_err_ene);
fprintf('The normalized error in determinat is %s.\n', norm_err_det);
fprintf('The normalized error in energy derivative (first variable) is %s.\n', norm_err_ene_d);
fprintf('The normalized error in determinant derivative (first variable) is %s.\n', norm_err_det_d);

%%%norm_err_ene in the order of 1e-8
%%%norm_err_det in the order of 1e-6
%%%norm_err_ened in the order of 1e-4
%%%norm_err_detd in the order of 1e-7


%%%Remark 1: Note that the real part of xdet is used for normalized error study. The complex number stems from
%%%SMW computations; see line 75 SMW_lkl. For better accuracy, the real part is considered 
%%%only at the end once all recursive computations are finished. So inside the code is
%%%not altered.

%%%Remark 2: In case you need derivative of the solution with respect to
%%%hyperparameters consider the following linear solve
%%%sol_d=-A^{-1}*(A_d_sol); This can be done by first using lkl_eval as
%%%above to find A_d_sol and then using A_d_sol as right hand side in one
%%%HMAT solve. HMAT solver is provided separately.




