function x = getX(par,utility,state,knext,munext,nunext,tauKbarnext,tauLbarnext,bnext,Approximator,g_grd,Pg)

beta  = par.beta;
delta = par.delta;
alpha = par.alpha;
Z     = par.z;

du    = utility.du;
duinv = utility.duinv;
ddu   = utility.ddu;

chi_tauL   = utility.chi_tauL;

dv    = @(l) utility.dv(l);

k=state.k;
g=state.g;
mu=state.mu;
tauLbar=state.tauLbar;
b=state.b;
nu=state.nu;

Eoutput=zeros(9,1);

for g_iter_next=1:length(g_grd)
    X=[knext;tauKbarnext;g_grd(g_iter_next);munext;tauLbarnext;nunext;bnext];
    Prediction=SimNetwork(X,Approximator.ymaxin,Approximator.yminin,Approximator.xmaxin,Approximator.xminin,Approximator.ymaxout,Approximator.yminout,Approximator.xmaxout,Approximator.xminout,Approximator.IW,Approximator.B1,Approximator.LW,Approximator.b2);
    Eoutput=Eoutput+Pg(state.g_iter,g_iter_next)*Prediction;
end


%Analytical stuff
c=duinv(beta*Eoutput(1)); % Investment Euler (mu)
l=((c+knext-(1-delta)*k+g)/(Z*k^alpha))^(1/(1-alpha)); % Resource constraint
Fk=Z*alpha*k^(alpha-1)*l^(1-alpha);
Fl=Z*k^alpha*(1-alpha)*l^(-alpha);
tauL=1-dv(l)/(Fl*du(c)); % Labor Supply (xi)
tauK=(b-beta*bnext*Eoutput(6)/du(c)+g-tauL*Fl*l)/((Fk-delta)*k); % Budget Constraint (nu)
xi=chi_tauL(tauL,tauLbar)/(Fl*du(c))-nunext*l; % FOC (tauL)
lambda=du(c)-munext*ddu(c)+mu*ddu(c)*(1+(Fk-delta)*(1-tauK))+nunext*ddu(c)*(tauK*(Fk-delta)*k+tauL*Fl*l-g-b)+nu*ddu(c)*b-xi*Fl*(1-tauL)*ddu(c); % FOC (c)

x(1)=c;
x(2)=l;
x(3)=knext;
x(4)=tauK;
x(5)=tauL;
x(6)=lambda;
x(7)=munext;
x(8)=nunext;
x(9)=xi;
x(10)=bnext;
x(11)=tauKbarnext;
x(12)=tauLbarnext;

end
