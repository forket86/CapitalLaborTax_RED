function residuals = FOCs5(par,utility,state,tauL,tauK,sigmanext,tauLbarnext,tauKbarnext,Approximator,g_grd,Pg,tauKmax,penalty)

beta  = par.beta;
delta = par.delta;
alpha = par.alpha;
gamma = par.gamma;
eta   = par.eta;
D     = par.D;
z     = par.z;
chi0  = par.chi0;
psi0  = par.psi0;

du    = utility.du;
ddu   = utility.ddu;


chi_tauL   = utility.chi_tauL;
psi_tauK   = utility.psi_tauK;

dv    = @(l) -utility.dv(l);
ddv   = @(l) -utility.ddv(l);

k=state.k;
tauKbar=state.tauKbar;
g=state.g;
tauLbar=state.tauLbar;

[l,c,knext,ltauL,ctauL,l_tauK,c_tauK,~,~] = getLaborCKnext(k,tauK,tauL,z,g,alpha,delta,gamma,eta,D);

Eoutput=zeros(5,1);
EDtauKbar=0;
EDtauLbar=0;
EDk=0;
Delta=1e-6;


for g_iter_next=1:length(g_grd)
    Eoutput=Eoutput+Pg(state.g_iter,g_iter_next)*SimNetwork([knext;tauKbarnext;g_grd(g_iter_next);tauLbarnext],Approximator.ymaxin,Approximator.yminin,Approximator.xmaxin,Approximator.xminin,Approximator.ymaxout,Approximator.yminout,Approximator.xmaxout,Approximator.xminout,Approximator.IW,Approximator.B1,Approximator.LW,Approximator.b2);
    
    H=Pg(state.g_iter,g_iter_next)*SimNetwork([knext+Delta;tauKbarnext;g_grd(g_iter_next);tauLbarnext],Approximator.ymaxin,Approximator.yminin,Approximator.xmaxin,Approximator.xminin,Approximator.ymaxout,Approximator.yminout,Approximator.xmaxout,Approximator.xminout,Approximator.IW,Approximator.B1,Approximator.LW,Approximator.b2);
    L=Pg(state.g_iter,g_iter_next)*SimNetwork([knext-Delta;tauKbarnext;g_grd(g_iter_next);tauLbarnext],Approximator.ymaxin,Approximator.yminin,Approximator.xmaxin,Approximator.xminin,Approximator.ymaxout,Approximator.yminout,Approximator.xmaxout,Approximator.xminout,Approximator.IW,Approximator.B1,Approximator.LW,Approximator.b2);
    EDk=EDk+(H(1)-L(1))/(2*Delta);
    
    H=Pg(state.g_iter,g_iter_next)*SimNetwork([knext;tauKbarnext+Delta;g_grd(g_iter_next);tauLbarnext],Approximator.ymaxin,Approximator.yminin,Approximator.xmaxin,Approximator.xminin,Approximator.ymaxout,Approximator.yminout,Approximator.xmaxout,Approximator.xminout,Approximator.IW,Approximator.B1,Approximator.LW,Approximator.b2);
    L=Pg(state.g_iter,g_iter_next)*SimNetwork([knext;tauKbarnext-Delta;g_grd(g_iter_next);tauLbarnext],Approximator.ymaxin,Approximator.yminin,Approximator.xmaxin,Approximator.xminin,Approximator.ymaxout,Approximator.yminout,Approximator.xmaxout,Approximator.xminout,Approximator.IW,Approximator.B1,Approximator.LW,Approximator.b2);
    EDtauKbar=EDtauKbar+(H(1)-L(1))/(2*Delta);
    
    H=Pg(state.g_iter,g_iter_next)*SimNetwork([knext;tauKbarnext;g_grd(g_iter_next);tauLbarnext+Delta],Approximator.ymaxin,Approximator.yminin,Approximator.xmaxin,Approximator.xminin,Approximator.ymaxout,Approximator.yminout,Approximator.xmaxout,Approximator.xminout,Approximator.IW,Approximator.B1,Approximator.LW,Approximator.b2);
    L=Pg(state.g_iter,g_iter_next)*SimNetwork([knext;tauKbarnext;g_grd(g_iter_next);tauLbarnext-Delta],Approximator.ymaxin,Approximator.yminin,Approximator.xmaxin,Approximator.xminin,Approximator.ymaxout,Approximator.yminout,Approximator.xmaxout,Approximator.xminout,Approximator.IW,Approximator.B1,Approximator.LW,Approximator.b2);
    EDtauLbar=EDtauLbar+(H(1)-L(1))/(2*Delta);
    
end

Eeuler1=Eoutput(1);
Eeuler2=Eoutput(2)-Eoutput(3);
EtauLnext=Eoutput(4);
EtauKnext=Eoutput(5);

lambda=(du(c)-sigmanext.*ddu(c).*knext-chi_tauL(tauL,tauLbar)./ctauL+((dv(l))).*ltauL./ctauL)./(1-(1-alpha).*z.*k.^(alpha).*l.^(-alpha).*ltauL./ctauL);
nu = du(c)-sigmanext.*ddu(c).*knext-lambda;
mu=(nu.*ctauL-chi_tauL(tauL,tauLbar))./ltauL;
%check1=-dv(l)-lambda*(1-alpha).*z.*k.^(alpha).*l.^(-alpha)-mu;
%check2=du(c)-sigmanext*ddu(c)*knext-nu-lambda;

%Implementability
residuals(1)=du(c).*knext-beta*Eeuler1;
%Eq. 17
residuals(2)=lambda-beta*Eeuler2+du(c).*sigmanext-sigmanext*beta*EDk;
%Eq. 18
residuals(3)=mu.*l_tauK-nu.*c_tauK+psi_tauK(tauK,tauKbar);
%residuals(7)=mu.*ltauL-nu.*ctauL+chi_tauL(tauL,tauLbar);
%residuals(8)=du(c)-sigmanext.*ddu(c).*knext-lambda-nu;
%Eq.23
residuals(4)=chi0*(EtauLnext-tauLbarnext)+sigmanext*EDtauLbar;
%Eq.24
residuals(5)=psi0*(EtauKnext-tauKbarnext)+sigmanext*EDtauKbar;