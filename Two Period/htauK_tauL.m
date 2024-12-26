function fp = htauK_tauL(tauL,k,g,z,alpha,D,eta)
%derivate of htauK function wrt tauL

%first solve l given k
l = hl(tauL,k,z,alpha,D,eta);
hlp = hl_tauL(tauL,k,z,alpha,D,eta);

fp = -(1-alpha)*hlp*g/(alpha*z*k^alpha*l^(2-alpha)) - (1-alpha)/alpha;
