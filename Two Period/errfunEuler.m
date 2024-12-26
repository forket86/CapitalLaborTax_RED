function [err,out] = errfunEuler(k,tauL,g,Pg,z,alpha,D,eta,beta)
%This function computes the error euler in the FC solution to the model, 
% for use in the solver find the optimal k for a given tauL choice

Ng = length(g);

%solve for l and tauK state by state
l = zeros(Ng,1);
tauK = zeros(Ng,1);
for i = 1:Ng
    l(i) = hl(tauL(i),k,g(i),z,alpha,D,eta);
    tauK(i) = htauK(tauL(i),k,g(i),z,alpha,D,eta);
end
y = z*k^alpha*l.^(1-alpha);
c = y - g;

%Euler expectation
E = alpha*z*k^(alpha-1)*(1-tauK).*l.^(1-alpha)./c;
E = E'*Pg;

%error in Euler
err = beta*E - 1;

%out
out.tauK = tauK;
out.l = l;
out.RHS = beta*E*k;