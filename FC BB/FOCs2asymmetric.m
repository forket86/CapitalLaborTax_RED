function residuals = FOCs2asymmetric(par,utility,state,tauL,tauK,sigmanext,tauLbarnext,tauKbarnext,Approximator,g_grd,Pg)

beta  = par.beta;
delta = par.delta;
alpha = par.alpha;
gamma = par.gamma;
eta   = par.eta;
D     = par.D;
z     = par.z;

du    = utility.du;
ddu   = utility.ddu;

chi_tauL   = utility.chi_tauL;
psi_tauK   = utility.psi_tauK;

dv    = @(l) -utility.dv(l);
ddv   = @(l) -utility.ddv(l);

k=state.k;
tauKbar=state.tauKbar;
g=state.g;
sigma=state.sigma;
tauLbar=state.tauLbar;

[l,c,knext,ltauL,ctauL,l_tauK,c_tauK,~,~] = getLaborCKnext(k,tauK,tauL,z,g,alpha,delta,gamma,eta,D);

Eoutput=zeros(5,1);
Echi_tauL=0;
Epsi_tauK=0;
for g_iter_next=1:length(g_grd)
    
    CurrentOutput=SimNetwork([knext;tauKbarnext;g_grd(g_iter_next);sigmanext;tauLbarnext],Approximator.ymaxin,Approximator.yminin,Approximator.xmaxin,Approximator.xminin,Approximator.ymaxout,Approximator.yminout,Approximator.xmaxout,Approximator.xminout,Approximator.IW,Approximator.B1,Approximator.LW,Approximator.b2);
    Eoutput=Eoutput+Pg(state.g_iter,g_iter_next)*CurrentOutput;
    
    Echi_tauL = Echi_tauL+Pg(state.g_iter,g_iter_next)*chi_tauL(CurrentOutput(4),tauLbarnext);
    Epsi_tauK = Epsi_tauK+Pg(state.g_iter,g_iter_next)*psi_tauK(CurrentOutput(5),tauKbarnext);
end


Eeuler1=Eoutput(1);
Eeuler2=Eoutput(2)-Eoutput(3);

lambda=(du(c)-(sigmanext-sigma).*ddu(c).*knext+sigma.*(du(c)+ddu(c).*c)-chi_tauL(tauL,tauLbar)./ctauL+((dv(l)+sigma.*(dv(l)+ddv(l).*l))).*ltauL./ctauL)./(1-(1-alpha).*z.*k.^(alpha).*l.^(-alpha).*ltauL./ctauL);
nu = du(c)-(sigmanext-sigma).*ddu(c).*knext+sigma.*(du(c)+ddu(c).*c)-lambda;
mu=(nu.*ctauL-chi_tauL(tauL,tauLbar))./ltauL;

%Implementability
residuals(1)=du(c).*knext-beta*Eeuler1;
%Eq. 15
residuals(2)=lambda-beta*Eeuler2+du(c).*(sigmanext-sigma);
%Eq. 16
residuals(3)=mu.*l_tauK-nu.*c_tauK+psi_tauK(tauK,tauKbar);
%Eq.18
residuals(4)=Echi_tauL;%EtauLnext-tauLbarnext;
%Eq.24
residuals(5)=Epsi_tauK;%EtauKnext-tauKbarnext;
%
% phi=10;
% if sigmanext<Sigma_L
%     residuals(6) = -phi*(Sigma_L-sigmanext)+log(1+phi*(Sigma_L-sigmanext));
% elseif sigmanext>Sigma_H
%     residuals(6) = phi*(sigmanext-Sigma_H)+log(1+phi*(sigmanext-Sigma_H));
% else
%     residuals(6) = 0;
% end