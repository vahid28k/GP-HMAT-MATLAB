
function [gp_H_mean,gp_H_var]=gp_H_eval(nodes_train,y,nodes_test,params,nmode)
warning('off');
%Case 1 reduces the kernel matrix, kl_tstr i.e. the kernel between the test and train data 
%Case 2 uses the full kernel between test and train data
delta1=params{4};
if params{3}>=size(nodes_test,2)
    nmode=2;
end

switch nmode
    case 1
        kernel_params=params{1};k=params{3};
        sol_y=back_solve(nodes_train,{y},params);
        [U,S,V]=randomized_svd_ID(nodes_test,nodes_train,kernel_params,k);
        gp_H_mean=(U*diag(S))*(V'*sol_y);
        
        [sol_K_trts]=back_solve(nodes_train,{V},params);
        right_mult=diag(S)*(V'*sol_K_trts)*(diag(S)'*U');
        gp_H_var=1+delta1-sum(U.*right_mult',2);
        
        
    case 2
        %K_tstr means test nodes are in the first dimension of matrix and
        %train nodes are in the second dimension
        
        kernel_params=params{1};
        sol_y=back_solve(nodes_train,{y},params);
        K_tstr=kernel_computation(nodes_test,nodes_train,kernel_params);
        gp_H_mean=K_tstr*sol_y;
        
        [sol_K_trts]=back_solve(nodes_train,{K_tstr'},params);
        gp_H_var=1+delta1-sum(K_tstr.*sol_K_trts',2);
end

