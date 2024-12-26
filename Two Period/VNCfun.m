function [VNC,out] = VNCfun(tauKbar,tauLbar,g,Pg,z,alpha,D,eta,beta,gammaK,gammaL,kFB)
%This function computes the value to a LTC planner of a certain choice of
%tauKbar. This is given to a minimiser to find the optimal LTC policy.

%solve Euler equation for optimal k
f = @(k) errfunEulerNC(k,tauKbar,tauLbar,g,Pg,z,alpha,D,eta,beta,gammaK,gammaL);
opts = optimset('TolFun',1e-10,'TolX',1e-10,'Display','off');
init = kFB;
[k,~,exf] = fsolve(f,init,opts);
if exf<1, error('solver'), end

%compute other variables
[~,out] = errfunEulerNC(k,tauKbar,tauLbar,g,Pg,z,alpha,D,eta,beta,gammaK,gammaL);
l = out.l;
tauK = out.tauK;
tauL = out.tauL;

%expectation in welfare
E = log(z*k^alpha*l.^(1-alpha) - g) - D*l.^(1+eta)/(1+eta) - gammaK/2*(tauK - tauKbar).^2 - gammaL/2*(tauL - tauLbar).^2;
E = E'*Pg;

%compute welfare
VNC = - k + beta*E;

%welfare excluding cost
VtilNC = VNC + (gammaK/2*(tauK' - tauKbar).^2 + gammaL/2*(tauL' - tauLbar).^2)*Pg;
out.VtilNC = VtilNC;
out.k = k;