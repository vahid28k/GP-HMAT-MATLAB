function [xdet,xdet_d,xene,xene_d,sol_out,A_d_sol]=SMW_lkl(indices_agg,offD_USV,offD_USVd,offD_dUSV,nodes,y,params)

floor_cut=params{6};

ncell=size(y,2);
ylength=zeros(1,ncell);
for j=1:ncell
    yj=y{1,j};
    ylength(j)=size(yj,2);
end
ranki=length(offD_USV{2});



n1=length(indices_agg{1});n2=length(indices_agg{2});
[n1cell,inds1]=inds_ncell(n1,floor_cut);y1=y(1:n1cell,:);
[n2cell,inds2]=inds_ncell(n2,floor_cut);y2=y(n1cell+1:n1cell+n2cell,:);


params_inds1={ylength,ranki,inds1};
params_inds2={ylength,ranki,inds2};


offD_U1=mat2cell(offD_USV{1},inds1);
offD_U2=mat2cell(offD_USV{3},inds2);



[xdet1,xdet1_d,xD1,ql1,qlr1,qry1,A_d_xD1,A_d_ql1,qlr_d1]=SMW_lkl_ingredients(nodes(:,indices_agg{1}),offD_U1,y1,params,params_inds1);
[xdet2,xdet2_d,xD2,ql2,qlr2,qry2,A_d_xD2,A_d_ql2,qlr_d2]=SMW_lkl_ingredients(nodes(:,indices_agg{2}),offD_U2,y2,params,params_inds2);



delta2=params{5};
corr=[diag(1./offD_USV{2}) qlr2;qlr1 diag(1./offD_USV{2})']+delta2*eye(2*ranki);
%corr_inv=inv(corr);%rank(corr)

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



%%%preparation of matrix & vector variables for energy and determinant
%%%computation
xD1_mat=cell2mat(xD1(:,1));
xD2_mat=cell2mat(xD2(:,1));
ql1_mat=cell2mat(ql1);
ql2_mat=cell2mat(ql2);
y1_mat=cell2mat(y1(:,1));
y2_mat=cell2mat(y2(:,1));


%%%%energy and determinant computation
xene11=y1_mat'*xD1_mat;
xene12=y2_mat'*xD2_mat;
xene21=[y1_mat'*ql1_mat,y2_mat'*ql2_mat];
xene23=[ql2_mat'*y2_mat;ql1_mat'*y1_mat];
sol_ene21=(corr'\xene21')';
sol_ene23=corr\xene23;
xene=(xene11+xene12)-(xene21*sol_ene23);

xdet3_aux=[offD_USV{2}];
xdet3=2*sum(log(xdet3_aux));
xdet4=sum((log(eig(corr))));
xdet=xdet1+xdet2+xdet3+xdet4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Derivative Business%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%preparation of variables
Udcell=offD_USVd{1};Sdcell=offD_USVd{2};Vdcell=offD_USVd{3};
dUcell=offD_dUSV{1};dScell=offD_dUSV{2};dVcell=offD_dUSV{3};
n_der=size(Udcell,2);
xdet_d=cell(1,n_der);
xene_d=cell(1,n_der);


%%Computation of A_d_sol, which is the result of multiplication of kernel
%%derivative to the solution of linear system characterized in the form of
%%SMW calculation i.e. A_d *(diagonal solution - corrective solution)
%%%A_d is a matrix with low rank off diagonals 
A_d_sol=cell(1,n_der);
for j=1:n_der
    
    A_d_xD1j=A_d_xD1{j};A_d_ql1j=A_d_ql1{j};
    A_d_xD2j=A_d_xD2{j};A_d_ql2j=A_d_ql2{j};
    Ud=Udcell{j};Sd=diag(Sdcell{j});Vd=Vdcell{j};
    A_d_sols=cell(1,size(xD1,2));
    
    for i=1:size(xD1,2)
        
        a11=cell2mat(A_d_xD1j(:,i))+(Ud*Sd)*(Vd'*cell2mat(xD2(:,i)));
        a12=cell2mat(A_d_ql1j(:,1))*sol_qry{i}(1:ranki,:)+(Ud*Sd)*(Vd'*ql2_mat)*sol_qry{i}(ranki+1:end,:);
        a21=(Vd*Sd)*(Ud'*cell2mat(xD1(:,i)))+cell2mat(A_d_xD2j(:,i));
        a22=(Vd*Sd)*(Ud'*ql1_mat)*sol_qry{i}(1:ranki,:)+cell2mat(A_d_ql2j(:,1))*sol_qry{i}(ranki+1:end,:);
        
        A_d_sols(i)={[a11-a12;a21-a22]};
                
    end
    
    A_d_sol(j)={A_d_sols};
    
end



for i=1:n_der
    A_d_xD1i=A_d_xD1{i};
    A_d_xD2i=A_d_xD2{i};
    %%%computation of intermediate quantities for derivative
    qlr1_d=2*dUcell{i}'*ql1_mat-qlr_d1{i};
    qlr2_d=2*dVcell{i}'*ql2_mat-qlr_d2{i};
    corr_d_aux=-dScell{i}./(offD_USV{2}.^2);
    corr_d=[diag(corr_d_aux) qlr2_d;qlr1_d diag(corr_d_aux)'];
    %%%Determinant derivative
    xdet3_d=2*sum(dScell{i}./xdet3_aux);
    xdet4_d=trace(corr\corr_d);
    xdet_d(i)={xdet1_d{i}+xdet2_d{i}+xdet3_d+xdet4_d};

    %%%%Energy derivative
    xene_d_aux1=xD1_mat'*dUcell{i}-cell2mat(A_d_xD1i(:,1))'*ql1_mat;
    xene_d_aux2=xD2_mat'*dVcell{i}-cell2mat(A_d_xD2i(:,1))'*ql2_mat;
    xene1_d = -xD1_mat'*cell2mat(A_d_xD1i(:,1))-xD2_mat'*cell2mat(A_d_xD2i(:,1));
    
    xene21_d = [xene_d_aux1,xene_d_aux2] ;
    xene23_d = [xene_d_aux2';xene_d_aux1'];
    xene22_d_augmented = -sol_ene21*corr_d*sol_ene23;
    xene_d(i) = {xene1_d - xene22_d_augmented - xene21_d*sol_ene23 - sol_ene21*xene23_d};
end



%%%%%%%%%%%%%%computation ingredients%%%%

function [xdet,xdet_d,xD,ql,qlr,qry,A_d_xD,A_d_ql,qlr_d]=SMW_lkl_ingredients(nodesi,offD_Ui,yi,params,params_inds)

kernel_params=params{1};
cutoff_size=params{2};
ylength=params_inds{1};
ranki=params_inds{2};
inds=params_inds{3};
delta1=params{4};
ni=size(nodesi,2);

if ni<=cutoff_size
    Aiicell=kernel_computation(nodesi,nodesi,kernel_params);
    Aii=Aiicell{1}+delta1*eye(ni);
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
    xdet=sum(log(eig(Aii)));
    
    
    n_der=size(Aiicell,2)-1;
    xdet_d=cell(1,n_der);A_d_xD=cell(1,n_der);A_d_ql=cell(1,n_der);
    qlr_d=cell(1,n_der);
    for j=2:size(Aiicell,2)
        A_d_sol_xD_ql=Aiicell{j}*sol_xD_ql;
        A_d_sol_xD_ql_cell=mat2cell(A_d_sol_xD_ql,size(A_d_sol_xD_ql,1),[ylength ranki]);
        A_d_xD(j-1)={A_d_sol_xD_ql_cell(:,1:end-1)};
        A_d_ql(j-1)={A_d_sol_xD_ql_cell(:,end)};
        xdet_d(j-1)={trace(Aii\Aiicell{j})};
        qlr_d(j-1)={cell2mat(ql)'*cell2mat(A_d_sol_xD_ql_cell(:,end))};
    end
    
else
    
    ncell=size(yi,2);
    rhs=yi;
    rhs(:,end+1)=offD_Ui;
    
    [xdet,xdet_d,~,~,sol_xD_ql,A_d_sol_xD_ql]=lkl_eval(nodesi,rhs,params);
    xD=(sol_xD_ql(:,1:end-1));
    ql=(sol_xD_ql(:,end));
    
    n_der=size(A_d_sol_xD_ql,2);
    A_d_xD=cell(1,n_der);A_d_ql=cell(1,n_der);
    qlr_d=cell(1,n_der);
    for j=1:n_der
        A_d_sol_xD_qlj=A_d_sol_xD_ql{j};
        A_d_xD(j)={A_d_sol_xD_qlj(:,1:end-1)};
        A_d_ql(j)={A_d_sol_xD_qlj(:,end)};
        qlr_d(j)={cell2mat(ql)'*cell2mat(A_d_sol_xD_qlj(:,end))};
    end
    
    
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




