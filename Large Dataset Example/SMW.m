function [sol_out]=SMW(indices_agg,offD_U,offD_S_inv,nodes,y,params)

floor_cut=params{6};

ncell=size(y,2);
ylength=zeros(1,ncell);
for j=1:ncell
    yj=y{1,j};
    ylength(j)=size(yj,2);
end
ranki=length(offD_S_inv{1});



n1=length(indices_agg{1});n2=length(indices_agg{2});
[n1cell,inds1]=inds_ncell(n1,floor_cut);y1=y(1:n1cell,:);
[n2cell,inds2]=inds_ncell(n2,floor_cut);y2=y(n1cell+1:n1cell+n2cell,:);


params_inds1={ylength,ranki,inds1};
params_inds2={ylength,ranki,inds2};


offD_U1=mat2cell(offD_U{1},inds1);
offD_U2=mat2cell(offD_U{2},inds2);


[xD1,ql1,qlr1,qry1]=SMW_ingredients(nodes(:,indices_agg{1}),offD_U1,y1,params,params_inds1);
[xD2,ql2,qlr2,qry2]=SMW_ingredients(nodes(:,indices_agg{2}),offD_U2,y2,params,params_inds2);

delta2=params{5};
corr=[offD_S_inv{1} qlr2;qlr1 offD_S_inv{1}']+delta2*eye(2*ranki);
%corr_inv=inv(corr);


sol_out=cell(size(ql1,1)+size(ql2,1),size(xD1,2));
sol_qry=cell(1,size(xD1,2));
for i=1:size(xD1,2)
    sol_qry(i)={corr\[qry2{i};qry1{i}]};
    for jj=1:size(ql1,1)
        %sol_out(jj,i)={xD1{jj,i}-ql1{jj}*[corr_inv(1:ranki,1:ranki)*qry2{i}+corr_inv(1:ranki,ranki+1:end)*qry1{i}]};
        sol_out(jj,i)={xD1{jj,i}-ql1{jj}*sol_qry{i}(1:ranki,:)};
    end
    for jj=1:size(ql2,1)
        %sol_out(jj+size(ql1,1),i)={xD2{jj,i}-ql2{jj}*[corr_inv(ranki+1:end,1:ranki)*qry2{i}+corr_inv(ranki+1:end,ranki+1:end)*qry1{i}]};
        sol_out(jj+size(ql1,1),i)={xD2{jj,i}-ql2{jj}*sol_qry{i}(ranki+1:end,:)};
    end
end


function [xD,ql,qlr,qry]=SMW_ingredients(nodesi,offD_Ui,yi,params,params_inds)

kernel_params=params{1};
cutoff_size=params{2};
ylength=params_inds{1};
ranki=params_inds{2};
inds=params_inds{3};
delta1=params{4};
ni=size(nodesi,2);

if ni<=cutoff_size
    Aii=kernel_computation(nodesi,nodesi,kernel_params)+delta1*eye(ni);
    rhs=yi;
    rhs(:,end+1)=offD_Ui;
    sol_xD_ql=Aii\cell2mat(rhs);
   
    sol_xD_ql_cell=mat2cell(sol_xD_ql,inds,[ylength ranki]);
    xD=sol_xD_ql_cell(:,1:end-1);
    ql=sol_xD_ql_cell(:,end);
    
    
    qry_qlr=cell2mat(offD_Ui)'*sol_xD_ql;
    qry_qlr_cell=mat2cell(qry_qlr,ranki,[ylength ranki]);
    qry=qry_qlr_cell(:,1:end-1);
    qlr=qry_qlr_cell{:,end};
    
else
    
    ncell=size(yi,2);
    rhs=yi;
    rhs(:,end+1)=offD_Ui;
    
    sol_xD_ql=back_solve(nodesi,rhs,params);
    xD=(sol_xD_ql(:,1:end-1));
    ql=(sol_xD_ql(:,end));
   
    qlr=zeros(ranki,ranki);
    for jj=1:size(ql,1)
        qlr=qlr+offD_Ui{jj}'*ql{jj};
    end
    qry=cell(1,ncell);
    for i=1:ncell
        qryi=zeros(ranki,ylength(i));
        for jj=1:size(ql,1)
            qryi=qryi+offD_Ui{jj}'*xD{jj,i};
        end
        qry(i)={qryi};
    end
    
end

function [ncell,inds]=inds_ncell(n,floor_cut)

if mod(n,floor_cut)==0
    ncell=n/floor_cut;
    inds=floor_cut*ones(1,ncell);
else
    ncell=floor(n/floor_cut)+1;
    inds=[floor_cut*ones(1,ncell-1) mod(n,floor_cut)];
end



