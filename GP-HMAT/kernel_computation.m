function [Kmat]=kernel_computation(x1,x2,kernel_params)
if kernel_params{3}==0
    nkernel=kernel_params{1};l=kernel_params{2};
    [Kmat]=kernel_generate_f(x1,x2,nkernel,l);
elseif kernel_params{3}==1
    if length(kernel_params)==3
        nkernel=kernel_params{1};l=kernel_params{2};
        [Kmat]=kernel_generate_fd(x1,x2,nkernel,l);
    else
        nkernel=kernel_params{1};l=kernel_params{2};
        [Kmat]=kernel_generate_fdi(x1,x2,nkernel,l,kernel_params{4});
    end
    
end

function [Kl]=kernel_generate_f(x1,x2,n,l)
sizel=(size(l,1)==1);

switch sizel
    
    case 1
        switch n
            case 1
                x1_x2=(repmat(x1(:),1,size(x2,2))-repmat(x2,size(x1,2),1)).^2/(l^2);
                Kl=reshape(exp(-sqrt(sum(reshape(x1_x2,size(x1,1),[]),1))),size(x1,2),[]);
            case 2
                x1_x2=(repmat(x1(:),1,size(x2,2))-repmat(x2,size(x1,2),1)).^2/(2*l^2);
                Kl=reshape(exp(-sum(reshape(x1_x2,size(x1,1),[]),1)),size(x1,2),[]);
        end
    case 0
        switch n
            case 1
                x1_x2=[(repmat(x1(:),1,size(x2,2))-repmat(x2,size(x1,2),1)).^2]./(repmat(l.^2,size(x1,2),size(x2,2)));
                Kl=reshape(exp(-sqrt(sum(reshape(x1_x2,size(x1,1),[]),1))),size(x1,2),[]);
            case 2
                x1_x2=[(repmat(x1(:),1,size(x2,2))-repmat(x2,size(x1,2),1)).^2]./(repmat(2*l.^2,size(x1,2),size(x2,2)));
                Kl=reshape(exp(-sum(reshape(x1_x2,size(x1,1),[]),1)),size(x1,2),[]);
        end
end


function [Kl_out]=kernel_generate_fd(x1,x2,n,l)
sizel=(size(l,1)==1);
switch sizel
    case 1
        switch n
            case 1
                x1_x2_open=(repmat(x1(:),1,size(x2,2))-repmat(x2,size(x1,2),1)).^2/(l^2);
                x1_x2=reshape(sqrt(sum(reshape(x1_x2_open,size(x1,1),[]),1)),size(x1,2),[]);
                Kl=exp(-x1_x2);
                Kl_d=(x1_x2/l).*Kl;
                Kl_out={Kl, Kl_d};
            case 2
                x1_x2_open=(repmat(x1(:),1,size(x2,2))-repmat(x2,size(x1,2),1)).^2/(2*l^2);
                x1_x2=reshape(sum(reshape(x1_x2_open,size(x1,1),[]),1),size(x1,2),[]);
                Kl=exp(-x1_x2);
                Kl_d=(x1_x2/(0.5*l)).*Kl;
                Kl_out={Kl, Kl_d};
        end
    case 0
        switch n
            case 1
                Kl_out=cell(1,size(x1,1)+1);
                x1_x2_open=[(repmat(x1(:),1,size(x2,2))-repmat(x2,size(x1,2),1)).^2]./(repmat(l.^2,size(x1,2),size(x2,2)));
                x1_x2=reshape(sqrt(sum(reshape(x1_x2_open,size(x1,1),[]),1)),size(x1,2),[]);
                Kl=exp(-x1_x2);
                Kl_out(1)={Kl};
                x1_x2(x1_x2==0)=1;
                Kl_d_tmp=x1_x2_open./repmat(l,size(x1,2),size(x2,2));
                for i=1:size(x1,1)
                    Kl_out(1+i)={(Kl_d_tmp(i:size(x1,1):end,:)./x1_x2).*Kl};
                end
            case 2
                Kl_out=cell(1,size(x1,1)+1);
                x1_x2_open=[(repmat(x1(:),1,size(x2,2))-repmat(x2,size(x1,2),1)).^2]./(repmat(2*l.^2,size(x1,2),size(x2,2)));
                x1_x2=reshape(sum(reshape(x1_x2_open,size(x1,1),[]),1),size(x1,2),[]);
                Kl=exp(-x1_x2);
                Kl_out(1)={Kl};
                Kl_d_tmp=x1_x2_open./repmat(0.5*l,size(x1,2),size(x2,2));
                for i=1:size(x1,1)
                    Kl_out(1+i)={Kl_d_tmp(i:size(x1,1):end,:).*Kl};
                end
        end
end


function [Kl_out]=kernel_generate_fdi(x1,x2,n,l,i)
sizel=(size(l,1)==1);
d=size(x1,1);
switch sizel
    case 1
        switch n
            case 1
                
                x1_x2_open=(repmat(x1(:),1,size(x2,2))-repmat(x2,size(x1,2),1)).^2/(l^2);
                x1_x2=reshape(sqrt(sum(reshape(x1_x2_open,size(x1,1),[]),1)),size(x1,2),[]);
                Kl=exp(-x1_x2);
                Kl_d=(x1_x2/l).*Kl;
                Kl_out={Kl, Kl_d};
            case 2
                x1_x2_open=(repmat(x1(:),1,size(x2,2))-repmat(x2,size(x1,2),1)).^2/(2*l^2);
                x1_x2=reshape(sum(reshape(x1_x2_open,size(x1,1),[]),1),size(x1,2),[]);
                Kl=exp(-x1_x2);
                Kl_d=(x1_x2/(0.5*l)).*Kl;
                Kl_out={Kl, Kl_d};
        end
    case 0
        switch n
            case 1
                Kl_out=cell(1,2);
                x1_x2_open=[(repmat(x1(:),1,size(x2,2))-repmat(x2,size(x1,2),1)).^2]./(repmat(l.^2,size(x1,2),size(x2,2)));
                x1_x2=reshape(sqrt(sum(reshape(x1_x2_open,size(x1,1),[]),1)),size(x1,2),[]);
                Kl=exp(-x1_x2);
                Kl_out(1)={Kl};
                x1_x2(x1_x2==0)=1;
                Kl_d_tmp=x1_x2_open(i:size(x1,1):end,:)./repmat(l(i,1),size(x1,2),size(x2,2));
                Kl_out(2)={(Kl_d_tmp./x1_x2).*Kl};
            case 2
                Kl_out=cell(1,2);
                x1_x2_open=[(repmat(x1(:),1,size(x2,2))-repmat(x2,size(x1,2),1)).^2]./(repmat(2*l.^2,size(x1,2),size(x2,2)));
                x1_x2=reshape(sum(reshape(x1_x2_open,size(x1,1),[]),1),size(x1,2),[]);
                Kl=exp(-x1_x2);
                Kl_out(1)={Kl};
                Kl_d_tmp=x1_x2_open(i:size(x1,1):end,:)./repmat(0.5*l(i,1),size(x1,2),size(x2,2));
                Kl_out(2)={Kl_d_tmp.*Kl};                 
        end
        
end


