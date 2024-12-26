function residuals = FOCs2(par,utility,state,knext,munext,nunext,tauKbarnext,tauLbarnext,bnext,Approximator,g_grd,Pg,phi,L,H)

beta  = par.beta;
delta = par.delta;
alpha = par.alpha;
Z     = par.z;

du    = utility.du;
duinv = utility.duinv;
ddu   = utility.ddu;

chi_tauL   = utility.chi_tauL;
psi_tauK   = utility.psi_tauK;

dv    = @(l) utility.dv(l);
ddv   = @(l) utility.ddv(l);

k=state.k;
tauKbar=state.tauKbar;
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

residuals=zeros(6,1);

%Analytical stuff
c=duinv(beta*Eoutput(1)); % Investment Euler (mu)
l=((c+knext-(1-delta)*k+g)/(Z*k^alpha))^(1/(1-alpha)); % Resource constraint
Fk=Z*alpha*k^(alpha-1)*l^(1-alpha);
Fl=Z*k^alpha*(1-alpha)*l^(-alpha);
Fll=Z*k^alpha*(1-alpha)*(-alpha)*l^(-alpha-1);
Fkl=Z*alpha*k^(alpha-1)*(1-alpha)*l^(-alpha);
tauL=1-dv(l)/(Fl*du(c)); % Labor Supply (xi)
tauK=(b-beta*bnext*Eoutput(6)/du(c)+g-tauL*Fl*l)/((Fk-delta)*k); % Budget Constraint (nu)
%nunexttemp=(mu+psi_tauK(tauK,tauKbar)/(du(c)*(Fk-delta)))/k; % FOC (tauK)
xi=chi_tauL(tauL,tauLbar)/(Fl*du(c))-nunext*l; % FOC (tauL)
lambda=du(c)-munext*ddu(c)+mu*ddu(c)*(1+(Fk-delta)*(1-tauK))+nunext*ddu(c)*(tauK*(Fk-delta)*k+tauL*Fl*l-g-b)+nu*ddu(c)*b-xi*Fl*(1-tauL)*ddu(c); % FOC (c)

% Resource constraint
%residuals(1)=c+knext-(1-delta)*k+g-F;

% Investment Euler (mu)
%residuals(2)=du(c)-beta*Eoutput(1);

% Budget Constraint (nu)
%residuals(3)=du(c)*(tauK*(Fk-delta)*k+tauL*Fl*l-g)+beta*bnext*Eoutput(6)-du(c)*b;

% Labor Supply (xi)
%residuals(4)=dv(l)-Fl*du(c)*(1-tauL);

% FOC (c)
%residuals(5)=du(c)-munext*ddu(c)+mu*ddu(c)*(1+(Fk-delta)*(1-tauK))+nunext*ddu(c)*(tauK*(Fk-delta)*k+tauL*Fl*l-g-b)+nu*ddu(c)*b-xi*Fl*(1-tauL)*ddu(c)-lambda;

% FOC (l)
residuals(1)=dv(l)-mu*du(c)*Fkl*(1-tauK)-nunext*du(c)*(tauK*Fkl*k+tauL*(Fll*l+Fl))-xi*(ddv(l)-Fll*du(c)*(1-tauL))-lambda*Fl;

% FOC (k)
residuals(2)=beta*Eoutput(2)+mu*beta*Eoutput(3)+beta*Eoutput(4)-beta*Eoutput(5)-lambda;

% FOC (tauK)
residuals(3)=-mu*du(c)*(Fk-delta)+nunext*du(c)*(Fk-delta)*k-psi_tauK(tauK,tauKbar);

% FOC (tauL)
%residuals(9)=nunext*du(c)*Fl*l+xi*Fl*du(c)-chi_tauL(tauL,tauLbar);

% FOC (tauKbar)
residuals(4)=tauKbarnext-Eoutput(9);

% FOC (tauLbar)
residuals(5)=tauLbarnext-Eoutput(8);

% FOC (b)
CN_up = 0; CN_down = 0;
if bnext<L
    CN_down = phi*(L-bnext)+log(1+phi*(L-bnext));
elseif bnext>H
    CN_up = phi*(bnext-H)+log(1+phi*(bnext-H));
end
residuals(6)=-nunext*Eoutput(6)+Eoutput(7)-CN_down+CN_up;

end
