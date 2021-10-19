function [offD_U,offD_S_inv]=compute_offD(nodes,indices_agg,kernel_params,k_min)

offD_U=cell(1,2);
offD_S_inv=cell(1,1);

k=min(k_min,min(length(indices_agg{1}),length(indices_agg{2})));
nmode=3;
[U,S,V]=offD_factorization(nodes(:,indices_agg{1}),nodes(:,indices_agg{2}),kernel_params,k,nmode);
offD_U(1,1)={U};offD_U(1,2)={V};
if nmode==1 || nmode==3
    offD_S_inv(1,1)={diag(1./S)};
elseif nmode==2
    offD_S_inv(1,1)={S};
end
    

%By default case 3 is used; see line 7. 
function [U,S,V]=offD_factorization(x1,x2,kernel_params,k,nmode)
switch nmode
    
    case 1
        
        A=kernel_computation(x1,x2,kernel_params);
        [U,S,V]=svd(A);
        S=diag(S);
        tol=1e-12;
        S=S(S>tol);
        k_eff=min(k,length(S));
        U=U(:,1:k_eff);
        S=S(1:k_eff,1);
        V=V(1:k_eff,:);
                
    case 2
        
        
        max_size_innprod=10*k;
        if size(x2,2)>max_size_innprod
            inds_rand=randsample(size(x2,2),max_size_innprod)';
            x2_subset=x2(:,inds_rand);
        else
            x2_subset=x2;
        end
        omega=(randn(size(x2_subset,2),k));
        y=kernel_computation(x1,x2_subset,kernel_params);
        y=y*omega;
        
        [Q,~]=qr(y,0);
        [~,~,IND]=qr(Q','vector');
      
        x0l=[x1(:,IND(1:k))] ;
       
        x0r=x0l;
        U=kernel_computation(x1,x0l,kernel_params);
        S=kernel_computation(x0r,x0l,kernel_params);
        V=kernel_computation(x0r,x2,kernel_params)';
        
        
    case 3
        
        [U,S,V]=randomized_svd_ID(x1,x2,kernel_params,k);
             
end


