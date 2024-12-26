function [err,out] = errfunNC(tauL,k,tauKbar,tauLbar,g,z,alpha,D,eta,gammaK,gammaL)
%error in optimal tauL choice in the NC model for a given choice of tax
%promises. used to solve for optimal tauL for a given tax promise
%this is evaluated for one single g value

%first solve l given k
l = hl(tauL,k,g,z,alpha,D,eta);
hlp = hl_tauL(l,k,g,z,alpha,D,eta);
tauK = htauK(tauL,k,g,z,alpha,D,eta);
%tauKp = htauK_tauL(tauL,k,g,z,alpha,D,eta);

y = z*k^alpha*l^(1-alpha);
c = y - g;

nu = gammaK*(tauK-tauKbar)/(alpha*y);

err = hlp*( (1-alpha)*z*k^alpha*l^(-alpha)*(1/c + nu*(alpha*tauK+(1-alpha)*tauL)) - D*l^eta ) + nu*(1-alpha)*z*k^alpha*l^(1-alpha) - gammaL*(tauL - tauLbar);

out.l = l;
out.tauK = tauK;
out.y = y;
out.c = c;