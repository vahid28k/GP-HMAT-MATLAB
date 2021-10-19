function [U,S,V,Udcell,Sdcell,Vdcell,dUcell,dScell,dVcell]=randomized_svd_ID_lkl(x1,x2,kernel_params,k)



max_size_constant=1e6;
max_size_innprod=floor(max_size_constant/size(x1,2));
size_innprod=min(max(2*k,max_size_innprod),10*k);

if size(x2,2)>size_innprod
    inds_rand=randsample(size(x2,2),size_innprod)';
    x2_subset=x2(:,inds_rand);
else
    x2_subset=x2;
end


A12cell=kernel_computation(x1,x2_subset,kernel_params);
size_cell=size(A12cell,2);

A12=A12cell{1};
omega=(randn(size(x2_subset,2),k));
A12_omega=A12*omega;
[Q,~]=qr(A12_omega,0);
[~,RQ,IND]=qr(Q','vector');
[~,INDt]=sort(IND,'ascend');
R11=RQ(:,1:k);R12=RQ(:,k+1:end);
T=(R11)\R12;
Mat_IT=[eye(k) T];
X=Mat_IT(:,INDt)';


%%%A_jrow is a row subset of A 
kernel_params(3)={0};
A_jrow=kernel_computation(x1(:,IND(1:k)),x2,kernel_params);
[Wstart,Rstart]=qr(A_jrow',0);
Z=X*Rstart';
[U,S,V_]=svd(Z,0);
V=(Wstart*V_);S=diag(S);



Sdcell=cell(1,size_cell-1);
Udcell=cell(1,size_cell-1);
Vdcell=cell(1,size_cell-1);

for jj=2:size_cell
    
    
    A12=A12cell{jj};
    A12_omega=A12*omega;
    [Q,~]=qr(A12_omega,0);
    [~,RQ,IND]=qr(Q','vector');
    [~,INDt]=sort(IND,'ascend');
    R11=RQ(:,1:k);R12=RQ(:,k+1:end);
    T=(R11)\R12;
    Mat_IT=[eye(k) T];
    X=Mat_IT(:,INDt)';
    
    kernel_params(3)={1};kernel_params(4)={jj-1};
    A_jrowcell=kernel_computation(x1(:,IND(1:k)),x2,kernel_params);
    A_jrowjj=A_jrowcell{2};
    
    [Wstartjj,Rstartjj]=qr(A_jrowjj',0);
    Zjj=X*Rstartjj';
    [Ujj,Sjj,Vjj_]=svd(Zjj,0);
    Vjj=(Wstartjj*Vjj_);
    Sjj=diag(Sjj);
    
    
    Udcell(jj-1)={Ujj};Sdcell(jj-1)={Sjj};Vdcell(jj-1)={Vjj};
    
    
    
end

[dUcell,dScell,dVcell]=SVD_derivative(U,S,V,Udcell,Sdcell,Vdcell);


function [dUcell,dScell,dVcell]=SVD_derivative(U,S,V,Udcell,Sdcell,Vdcell)

%%This function takes the SVD of the kernel matrix and kernel derivative
%%matrix and produces the derivative of the kernel SVD. The analytical expressions are
%%from the paper below.


%https://j-towns.github.io/papers/svd-derivative.pdf
% [1] Alan Edelman, Tomas A Arias, and Steven T Smith. The geometry of algorithms
% with orthogonality constraints. SIAM journal on Matrix Analysis and Applications,
% 20(2):303-353, 1998.
% [2] Thomas P. Minka. Old and new matrix algebra useful for statistics.
% http://research.microsoft.com/en-us/um/people/minka/papers/matrix/minka-matrix.pdf.


%%%U,S,V is the SVD of the kernel
%%%Ud,Sd,Vd is the SVD of the kernel derivative
%%%dU,dS,dV is the derivative of the kernel SVD


S_orig=repmat(S',length(S),1).^2-repmat(S,1,length(S)).^2;
S_orig(S_orig==0)=1e20;
F=1./S_orig;F(abs(F)<=1e-18)=0;

size_cell=size(Udcell,2);

dScell=cell(1,size_cell);
dUcell=cell(1,size_cell);
dVcell=cell(1,size_cell);

for i=1:size_cell
    
    Ud=Udcell{i};Sd=Sdcell{i};Vd=Vdcell{i};
    
    UUd=U'*Ud;VVd=V'*Vd;
    USVd=UUd*diag(Sd)*VVd';
    Uaux=USVd*diag(S)+diag(S)*USVd';
    Vaux=diag(S)*USVd+USVd'*diag(S);
    
    dS=UUd*diag(Sd)*VVd';
    dU=[U*(F.*Uaux - USVd*diag(1./S) )] + [Ud*(diag(Sd)*VVd'*diag(1./S))];
    dV=[V*(F.*Vaux - USVd'*diag(1./S) )] + [Vd*(diag(Sd)*UUd'*diag(1./S))];
    
    dScell(i)={diag(dS)};dUcell(i)={dU};dVcell(i)={dV};
    
   
end




