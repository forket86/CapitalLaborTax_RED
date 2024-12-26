function fp = hl_tauL(l,k,g,z,alpha,D,eta)
%derivative of hl function wrt tauL

fp = -(1-alpha)*z*k^alpha/( (1+eta)*D*z*k^alpha*l^eta - (eta+alpha)*D*g*l^(eta+alpha-1) );
