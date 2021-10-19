clc;clear all;close all;

load data1e5.mat;

for i=1:4
    X(i,:)=X(i,:)/max(X(i,:));
    Y=Y/max(Y);
end

xtrain=X(:,1:9e4);ytrain=Y(1:9e4);m_y=mean(ytrain);
ytrain=ytrain-m_y;

xtest=X(:,9e4+1:end);ytest=Y(9e4+1:end);

global liter fiter cc;
liter=zeros(4,10000);fiter=zeros(1,10000);cc=1;

%%%In case you have vector valued y in d dimension 
%%%then ytrain will be given in the form of nsamples*d matrix.


options = optimoptions('fminunc','Display','iter','Algorithm','trust-region','SpecifyObjectiveGradient',true);
nkernel=2;delta1=1e-3;delta2=0e-8;n=size(xtrain,2);
k=20;cutoff_size=1005;
x0 = ones(4,1);
params={{nkernel,x0,1},cutoff_size,k,delta1,delta2};
fun = @(lp) gp_H_opt(xtrain,ytrain,params,lp);
tic;
lopt = fminunc(fun,x0,options);
t=toc;

%%In the case of vector valued y, the esitmate for variance does not depend
%%on the dimensionality of y, therefore the variances are scalar values at
%%each uknown node. 
params={{nkernel,lopt,0},cutoff_size,k,delta1,delta2};
[gp_mean,gp_var] = gp_H_eval(xtrain,ytrain,xtest,params,1);
gp_mean=gp_mean+m_y;


%k=20;
% Iteration        f(x)          step          optimality   CG-iterations
%      0            -226695                      1.24e+07                
%      1            -226695             10       1.24e+07           1
%      2            -227190            2.5       3.36e+07           0
%      3            -227198          0.625            108           1
%      4            -227198        0.15625            108           1
%      5            -227198      0.0390625            108           0
%      6            -227198     0.00976562            108           0
%      7            -227228     0.00244141       6.55e+04           0
%      8            -227228    0.000610352       6.55e+04           2
%      9            -227228    0.000152588       6.55e+04           0
%     10            -227228     3.8147e-05       6.55e+04           0
%     11            -227228    9.53674e-06       6.55e+04           0
%     12            -227228    2.38419e-06       6.55e+04           0
%     13            -227252    5.96046e-07       1.92e+04           0

%%lopt=[0.257671233575062;-1.14557764418883;0.260236864715620;-0.342555337673058];

%k=30
% Iteration        f(x)          step          optimality   CG-iterations
%      0            -227305                           384                
%      1            -227305             10            384           1
%      2            -227305            2.5            384           0
%      3            -227305          0.625            384           0
%      4            -227306        0.15625            430           0
%      5            -227311      0.0390625            575           1
%      6            -227311     0.00976562            575           2
%      7            -227311     0.00244141            575           0
%      8            -227311    0.000610352            575           0
%      9            -227311    0.000152588            575           0
%     10            -227311     3.8147e-05            575           0
%     11            -227311    9.53674e-06            575           0
%     12            -227311    2.38419e-06            575           0
%     13            -227311    5.96046e-07            575           0

%%lopt=[0.990709525357952;0.970458675753159;0.973239985880320;0.888236475685549];
