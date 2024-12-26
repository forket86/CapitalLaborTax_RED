function l = hl(tauL,k,g,z,alpha,D,eta)
%function solves for equilibrium labor supply given taxes and capital

%first solve l given k (same in all states in FB problem
f = @(x) (1-alpha)*z*k^alpha*(1-tauL).*exp(x).^(-alpha)./(z*k^alpha*exp(x).^(1-alpha) - g) - D*exp(x).^eta;

opts = optimset('TolFun',1e-10,'TolX',1e-10,'Display','off');
init = 0;
[x,~,exf] = fsolve(f,init,opts);
if exf<1, error('solver'), end
l = exp(x);
