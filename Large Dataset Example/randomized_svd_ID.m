function [U,S,V]=randomized_svd_ID(x1,x2,kernel_params,k)

%%%%max_size_constant can be changed to a higher number in case the number
%%%%of nodes are higher than 1e6. If the matrixis divided to 10 slices and
%%%%the number of nodes are 1e7 then every block will have 1e6 in the first
%%%%dimension so with max_size_constant=1e6 only one column will be
%%%%evaluated. This e.g. can change to 5e6 and and line 11 can change to size_innprod=min(max(5,max_size_innprod),10*k);
%%%%in that case at least 5 columns will be evaluated.


max_size_constant=1e6;
max_size_innprod=floor(max_size_constant/size(x1,2));
size_innprod=min(max(2*k,max_size_innprod),10*k);

if size(x2,2)>size_innprod
    inds_rand=randsample(size(x2,2),size_innprod)';
    x2_subset=x2(:,inds_rand);
else
    x2_subset=x2;
end


A12=kernel_computation(x1,x2_subset,kernel_params);
omega=(randn(size(x2_subset,2),k));
A12_omega=A12*omega;
[Q,~]=qr(A12_omega,0);


%Subspace Iteration
% qq=5;
% for i=1:qq
%     [Q,~]=qr(A12'*Q,0);
%     [Q,~]=qr(A12*Q,0);
% end



%%%%%ID%%%%%
[~,RQ,IND]=qr(Q','vector');
[~,INDt]=sort(IND,'ascend');

R11=RQ(:,1:k);R12=RQ(:,k+1:end);
T=(R11)\R12;
Mat_IT=[eye(k) T];
X=Mat_IT(:,INDt)';

%%%A_jrow is a row subset of A 
A_jrow=kernel_computation(x1(:,IND(1:k)),x2,kernel_params);
[Wstart,Rstart]=qr(A_jrow',0);
Z=X*Rstart';
[U,S,V_]=svd(Z,0);
V=(Wstart*V_);S=diag(S);









