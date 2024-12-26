function tauK = htauK(tauL,k,g,z,alpha,D,eta)
%function solves equilibrium tauK given tauL and k and state

%first solve l given k (same in all states in FB problem
l = hl(tauL,k,g,z,alpha,D,eta);

tauK = g/(alpha*z*k^alpha*l^(1-alpha)) - (1-alpha)/alpha*tauL;
