function [sol_out]=back_solve(nodes,y,params)
%parameters are: {nkernel l 0-1} cutoff_size k delta1 delta2

kernel_params=params{1};
cutoff_size=params{2};
k_uniform=params{3};


if size(y,2)==1
    
    [~,floor_cut]=matrix_slice_size(cutoff_size,k_uniform);
    floor_cut=min(100*floor(cutoff_size/100),floor_cut);
    params(6)={floor_cut};
    
    [ind_perm]=perm_generator(nodes,kernel_params,cutoff_size,[1:size(nodes,2)]);
    nodes=nodes(:,ind_perm);
    y=y{:,1}(ind_perm,:);
    [~,inds]=inds_ncell(size(y,1),floor_cut);
    y=mat2cell(y,inds);
    
end

floor_cut=params{6};
[k_min,slice_size]=matrix_slice_size(size(nodes,2),k_uniform);
indices_agg(1)={1:floor_cut*floor(slice_size/floor_cut)};
indices_agg(2)={floor_cut*floor(slice_size/floor_cut)+1:size(nodes,2)};



[offD_U,offD_S_inv]=compute_offD(nodes,indices_agg,kernel_params,k_min);

[sol_out]=SMW(indices_agg,offD_U,offD_S_inv,nodes,y,params);

if size(y,2)==1
    sol_out_mat=cell2mat(sol_out);
    sol_out_mat_perm(ind_perm,:)=sol_out_mat;
    sol_out=sol_out_mat_perm;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [k_min,slice_size]=matrix_slice_size(n,k_uniform)

if 1e5<n && n<=1e6
    slice_size=10^(floor(log10(n-0.5)));
    k_min=k_uniform;
elseif 1e4<n && n<=1e5
    slice_size=10^(floor(log10(n-0.5)));
    k_min=k_uniform;
elseif 1e3<n && n<=1e4
    slice_size=10^(floor(log10(n-0.5)));
    k_min=k_uniform;
elseif 1e2<n && n<=1e3
    slice_size=10^(floor(log10(n-0.5)));
    k_min=k_uniform;
elseif 1e1<n && n<=1e2
    slice_size=10^(floor(log10(n-0.5)));
    k_min=k_uniform;
end


function [ind_fin]=perm_generator(nodes,kernel_params,cutoff_size,ind_ini)

[~,slice_size]=matrix_slice_size(size(nodes,2),20);
[indices_agg]=permutation(nodes,kernel_params,slice_size);

ind_fin_1=perm_cutoff(nodes(:,indices_agg{1}),kernel_params,cutoff_size,ind_ini(indices_agg{1}));
ind_fin_2=perm_cutoff(nodes(:,indices_agg{2}),kernel_params,cutoff_size,ind_ini(indices_agg{2}));

ind_fin=[ind_fin_1 ind_fin_2];

function [indices_agg]=permutation(nodes,kernel_params,slice_size)
indices_agg=cell(1,2);
[Kl]=kernel_computation(nodes(:,1),nodes,kernel_params);
[~,p]=sort(Kl,'descend');
indices_agg(1)={p(1:slice_size)};
indices_agg(2)={p(slice_size+1:end)};

function [indices_fin]=perm_cutoff(nodes,kernel_params,cutoff_size,indices)
if size(nodes,2)<=cutoff_size
    indices_fin=indices;
else
    indices_fin=perm_generator(nodes,kernel_params,cutoff_size,indices);
end


function [ncell,inds]=inds_ncell(n,floor_cut)

if mod(n,floor_cut)==0
    ncell=n/floor_cut;
    inds=floor_cut*ones(1,ncell);
else
    ncell=floor(n/floor_cut)+1;
    inds=[floor_cut*ones(1,ncell-1) mod(n,floor_cut)];
end


