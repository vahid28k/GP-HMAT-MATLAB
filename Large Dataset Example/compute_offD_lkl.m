function [offD_USV,offD_USVd,offD_dUSV]=compute_offD_lkl(nodes,indices_agg,kernel_params,k_min)

k=min(k_min,min(length(indices_agg{1}),length(indices_agg{2})));
[U,S,V,Udcell,Sdcell,Vdcell,dUcell,dScell,dVcell]=randomized_svd_ID_lkl(nodes(:,indices_agg{1}),nodes(:,indices_agg{2}),kernel_params,k);

%%other options for compute off D
%[U,S,V,Udcell,Sdcell,Vdcell,dUcell,dScell,dVcell]=usual_svd(nodes(:,indices_agg{1}),nodes(:,indices_agg{2}),kernel_params,k);
%[U,S,V,Udcell,Sdcell,Vdcell,dUcell,dScell,dVcell]=randomized_svd_ID_lkl_alternative(nodes(:,indices_agg{1}),nodes(:,indices_agg{2}),kernel_params,k);


offD_USV={U,S,V};
offD_USVd={Udcell,Sdcell,Vdcell};
offD_dUSV={dUcell,dScell,dVcell};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Usual SVD%%%%%%%%%%%%%%%%%%%%%%%


function [U,S,V,Udcell,Sdcell,Vdcell,dUcell,dScell,dVcell]=usual_svd(x1,x2,kernel_params,k)
A12cell=kernel_computation(x1,x2,kernel_params);

A12=A12cell{1};
[U,S,V]=svd(A12);S=diag(S);S=S(1:k,1);U=U(:,1:k);V=V(:,1:k);


size_cell=size(A12cell,2);

Sdcell=cell(1,size_cell-1);
Udcell=cell(1,size_cell-1);
Vdcell=cell(1,size_cell-1);

for jj=2:size_cell
    
    
    A12=A12cell{jj};
    [Ujj,Sjj,Vjj]=svd(A12,0);
    Sjj=diag(Sjj);Sjj=Sjj(1:k,1);Ujj=Ujj(:,1:k);Vjj=Vjj(:,1:k);
    
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


