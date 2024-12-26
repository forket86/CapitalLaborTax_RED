%This file solves the two period model in the FC and LTC cases

clear all
close all
clc


%controls
load_init = 1; %load initial guesses or not

%parameters
alpha = 1/3;
eta = 2;
z = 1;
beta = 0.95;

%max and min tauL values
tauLmin = 0;
tauLmax = eta/(1+eta);

%g shocks probs
Ng = 2;
Pg = ones(Ng,1)/Ng;
gbarfrac = 0.075;
Deltag = 0.25;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build gamma grid and load initial values


%grid of gamma values -- the FC and LTC problems are solved for each of
%these values
gammaKgrid = [linspace(0,0.2,20),linspace(0.2,2,15)]';
gammaKgrid = unique(gammaKgrid);
Ngamma = length(gammaKgrid);

%imposed lower and upper bounds on capital tax promise in LTC solver (true 
%solution is interior to these values)
tauKbarmin = -1;
tauKbarmax = 0.05;

if load_init
    load initvals
else
    %FC with csc initial guess
    tauLFCc_init = 0.5*tauLmax*ones(Ngamma,Ng);
    %LTC with csc initial guess
    tauKbarNCc_init = (tauKbarmin+tauKbarmax)/2*ones(Ngamma,1);
end






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First solve FB problem

%solve for FB k and calibrate D so that l=1 in FB problem
f = @(x) errfunFB(exp(x(1:Ng)),exp(x(Ng+1)),exp(x(Ng+2)),z,alpha,eta,beta,Pg,gbarfrac,Deltag);
opts = optimset('TolFun',1e-10,'TolX',1e-10,'Display','off');
init = zeros(Ng+2,1);
[x,~,exf] = fsolve(f,init,opts);
if exf<1, error('solver'), end
lFB = exp(x(1:Ng));
kFB = exp(x(Ng+1));
D = exp(x(Ng+2));

%run again to get out l
[~,out] = errfunFB(lFB,kFB,D,z,alpha,eta,beta,Pg,gbarfrac,Deltag);
yFB = out.y;
cFB = out.c;
g = out.g;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve NC problem with no costly state contingency

%NC sets labour taxes to zero
tauLNC = zeros(Ng,1);
%solve kNC
f = @(k) errfunEuler(k,tauLNC,g,Pg,z,alpha,D,eta,beta);
opts = optimset('TolFun',1e-10,'TolX',1e-10,'Display','off');
init = kFB;
[kNC,~,exf] = fsolve(f,init,opts);
if exf<1, error('solver'), end

%compute other NC objects
[err,out] = errfunEuler(kNC,tauLNC,g,Pg,z,alpha,D,eta,beta);
lNC = out.l;
tauKNC = out.tauK;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve FC problem with no costly state contingency


%solve for no cost case first
gammaK = 0;
tauKbarFC = 0;
gammaL = 0;
tauLbarFC = 0;

%maximise value function
lb = zeros(Ng,1); %lower bound on labour tax
ub = tauLmax*ones(Ng,1); %upper bound
A = []; b = []; Aeq = []; beq = []; %redundant constraints
x0 = (lb + ub)/2; %initial guess
f = @(tauL) - VFCfun(tauL,g,Pg,z,alpha,D,eta,beta,gammaK,gammaL,kFB);
opts = optimset('TolFun',1e-10,'TolX',1e-10,'Display','off');
[tauLFC,fval,exf] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,[],opts);
if exf<1, error('solver'), end
%compute the rest
[VFC,out] = VFCfun(tauLFC,g,Pg,z,alpha,D,eta,beta,gammaK,gammaL,kFB);
tauKFC = out.tauK;
lFC = out.l;
kFC = out.k;
cFC = out.c;
yFC = out.y;
VtilFC = out.VtilFC;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve FC problem with costly state contingency


tauLFCc = zeros(Ngamma,Ng);
tauKFCc = zeros(Ngamma,Ng);
lFCc = zeros(Ngamma,Ng);
tauKbarFCc = zeros(Ngamma,1);
VFCc = zeros(Ngamma,1);
VtilFCc = zeros(Ngamma,1);
kFCc = zeros(Ngamma,1);




%loop over gamma values and solve the FC model for each value
for i = 1:Ngamma
    
    gammaK = gammaKgrid(i)
    gammaL = 0; %assuming no cost on labor tax
    
    %maximise value function
    lb = zeros(Ng,1); %lower bound on labour tax
    ub = tauLmax*ones(Ng,1); %upper bound
    A = []; b = []; Aeq = []; beq = []; %redundant constraints
    %x0 = (lb + ub)/2; %initial guess
    x0 = tauLFCc_init(i,:)';
    f = @(tauL) - VFCfun(tauL,g,Pg,z,alpha,D,eta,beta,gammaK,gammaL,kFB);
    opts = optimset('TolFun',1e-10,'TolX',1e-10,'Display','off');
    [x,fval,exf] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,[],opts);
    if exf<1, error('solver'), end
    %policies
    tauLFCc(i,:) = x;
    %compute the rest
    [VFCc_,out] = VFCfun(tauLFCc(i,:)',g,Pg,z,alpha,D,eta,beta,gammaK,gammaL,kFB);
    VFCc(i) = VFCc_;
    tauKFCc(i,:) = out.tauK;
    lFCc(i,:) = out.l';
    kFCc(i) = out.k;
    cFCc(i,:) = out.c';
    yFCc(i,:) = out.y';
    VtilFCc(i) = out.VtilFC;
    tauKbarFCc(i) = out.tauKbar;
    tauLbarFCc(i) = out.tauLbar;
    
end


%when gamma=0 just returns the initial value. true limit is well defined
%and so just skip this entry from graphs
tauKbarFCc(1) = NaN;

%save solution as initial guess
tauLFCc_init = tauLFCc;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve LTC (a.k.a NC) problem with costly state contingency

%solver parameters
lb = tauKbarmin; %lower bound on promise
ub = tauKbarmax; %upper bound
A = []; b = []; Aeq = []; beq = []; %redundant constraints

tauLNCc = zeros(Ngamma,Ng);
tauKNCc = zeros(Ngamma,Ng);
lNCc = zeros(Ngamma,Ng);
tauKbarNCc = zeros(Ngamma,1);
VNCc = zeros(Ngamma,1);
VtilNCc = zeros(Ngamma,1);
kNCc = zeros(Ngamma,1);


%loop over gamma values and solve the LTC model for each value
for i = 1:Ngamma
    
    gammaK = gammaKgrid(i)
    gammaL = 0;
    
    %maximise value function
    f = @(x) - VNCfun(x,0,g,Pg,z,alpha,D,eta,beta,gammaK,gammaL,kFB);
    opts = optimset('TolFun',1e-5,'TolX',1e-5,'Display','off');
    x0 = tauKbarNCc_init(i);
    [x,fval,exf] = fmincon(f,x0,A,b,Aeq,beq,lb,ub,[],opts);
    if exf<1, error('solver'), end
    tauKbarNCc(i) = x;
    %compute the rest
    [VNC_,out] = VNCfun(tauKbarNCc(i),0,g,Pg,z,alpha,D,eta,beta,gammaK,gammaL,kFB);
    VNCc(i) = VNC_;
    tauKNCc(i,:) = out.tauK;
    tauLNCc(i,:) = out.tauL;
    lNCc(i,:) = out.l';
    kNCc(i) = out.k;
    VtilNCc(i) = out.VtilNC;
    
    
end

%save initial value
tauKbarNCc_init = tauKbarNCc;

%when gamma=0 just returns the initial value. true limit is well defined
%and so just skip this entry from graphs
tauKbarNCc(1) = NaN;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save workspace

save results2P

%save initial guesses
save initvals tauLFCc_init tauKbarNCc_init



