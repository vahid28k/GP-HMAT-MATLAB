function [f,gvec] = gp_H_opt(nodes,y,params,l)
global liter fiter cc;
dy=size(y,2);

warning('off');
params{1}(2)={l};
n=size(nodes,2);
[xdet,xdet_d,energy,energy_d,~,~]=lkl_eval(nodes,{y},params);
f=-0.5*dy*xdet-0.5*trace(energy)-0.5*n*log(2*pi);
gvec=zeros(1,length(l));
for i=1:length(l)
    gvec(i)=-0.5*dy*xdet_d{i}-0.5*trace(energy_d{i});
end
f=-f;
gvec=-gvec;

liter(:,cc)=l;fiter(1,cc)=f;
cc=cc+1;