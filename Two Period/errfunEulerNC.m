function [err,out] = errfunEulerNC(k,tauKbar,tauLbar,g,Pg,z,alpha,D,eta,beta,gammaK,gammaL)
%error function to solve the Euler equation for optimal k value in the NC 
% model where the taxes are chosen ex-post as period 2 optimal values

Ng = length(g);

%solve for l and tauK state by state
l = zeros(Ng,1); c = zeros(Ng,1);
tauK = zeros(Ng,1);
tauL = zeros(Ng,1);
opts = optimset('TolFun',1e-10,'TolX',1e-10,'Display','off');
for i = 1:Ng
    f = @(x) errfunNC(exp(x),k,tauKbar,tauLbar,g(i),z,alpha,D,eta,gammaK,gammaL);
    init = -10;
    [x,~,exf] = fsolve(f,init,opts);
    if exf<1, error('solver'), end
    tauL(i) = exp(x);
    [~,out] = errfunNC(tauL(i),k,tauKbar,tauLbar,g(i),z,alpha,D,eta,gammaK,gammaL);
    l(i) = out.l;
    c(i) = out.c;
    tauK(i) = out.tauK;
end


%Euler expectation
E = alpha*z*k^(alpha-1)*(1-tauK).*l.^(1-alpha)./c;
E = E'*Pg;

%error in Euler
err = beta*E - 1;

%out
clear out
out.tauK = tauK;
out.tauL = tauL;
out.l = l;