function [err,out] = errfunFB(l,k,D,z,alpha,eta,beta,Pg,gbarfrac,Deltag)
%This function computes the error in the first best solution of the model,
%to be used in the solver to solve for the first best

%compute y given k,l guess
y = z*k^alpha*l.^(1-alpha);
%averages
Ey = y'*Pg;
El = l'*Pg;

%average g -- choose so that kNC = ~ 0.5*kFB. IMPORTANT: if g is set
%too high then no solution will exist in NC case since zero labour taxes
%will require too high capital taxes so no investment happens.
gbar = gbarfrac*Ey;
g = gbar*[1-Deltag;1+Deltag];

%c given y and g
c = y - g;

%first solve l given k (same in all states in FB problem
errl = (1-alpha)*z*k^alpha*l.^(-alpha)./c - D*l.^eta;

%now check error in Euler
errK = (beta*alpha*z*k^(alpha-1)*l.^(1-alpha)./c)'*Pg - 1;

err = [errK;
       errl;
       El-1];
   
out.g = g;
out.c = c;
out.y = y;