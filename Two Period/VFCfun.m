function [VFC,out] = VFCfun(tauL,g,Pg,z,alpha,D,eta,beta,gammaK,gammaL,kFB)
%This function computes the value to a FC planner of a certain choice of
%tauL. This is given to a minimiser to find the optimal FC policy.

%solve the Euler equation for optimal k
f = @(k) errfunEuler(k,tauL,g,Pg,z,alpha,D,eta,beta);
opts = optimset('TolFun',1e-10,'TolX',1e-10,'Display','off');
init = kFB;
[k,~,exf] = fsolve(f,init,opts);
if exf<1, error('solver'), end

%compute other variables
[~,out] = errfunEuler(k,tauL,g,Pg,z,alpha,D,eta,beta);
l = out.l;
tauK = out.tauK;

%expected policies and hence optimal announcements (FC so unbiased)
EtauK = tauK'*Pg;
EtauL = tauL'*Pg;
tauKbar = EtauK;
tauLbar = EtauL;

y = z*k^alpha*l.^(1-alpha);
c = y - g;

%expectation in welfare
E = log(c) - D*l.^(1+eta)/(1+eta) - gammaK/2*(tauK - tauKbar).^2 - gammaL/2*(tauL - tauLbar).^2;
E = E'*Pg;

%compute welfare
VFC = - k + beta*E;

%welfare excluding cost
VtilFC = VFC + (gammaK/2*(tauK' - tauKbar).^2 + gammaL/2*(tauL' - tauLbar).^2 )*Pg;
out.VtilFC = VtilFC;
out.k = k;
out.y = y;
out.c = c;
out.tauKbar = tauKbar;
out.tauLbar = tauLbar;