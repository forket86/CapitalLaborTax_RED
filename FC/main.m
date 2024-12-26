clear all;
clc;

addpath('../Utils/');

%% Parameters
beta=0.96;%discount factor
delta=0.08;%depreciation
alpha=0.36;%capital share
gamma=1.0;
eta=2.0;
T=10000; %time horizon
chi0=0.05;
psi0=100;
GPU=false;

% Stoch. params of g process
rho= 0.9774; %0.7785;
gnum=2;
sigmasq= (4/4)*0.0161; %0.0114; %0.07

rng(1)

%% Calibration to FC Steady State
% Targets
kFC = 1;
lFC = 1;
gfracFC = 0.2019;
%derived params
Z = (1/beta + delta - 1)/(alpha*kFC^(alpha-1)*lFC^(1-alpha) );
yFC = Z*kFC^alpha*lFC^(1-alpha);
g = gfracFC*yFC;
tauKFC=0;
tauLFC = g/((1-alpha)*yFC);
cFC = yFC - g - delta*kFC;
D = (1-tauLFC)*(1-alpha)*Z*kFC^alpha*lFC^(-alpha)/( lFC^eta*cFC^gamma );
z_grd=Z;

% Parameters structure
par.alpha = alpha;
par.beta = beta;
par.delta = delta;
par.sigma = gamma;
par.gamma = gamma;
par.eta = eta;
par.z = Z;
par.g = g;
par.D = D;
par.dd = 1;
par.chi0=chi0;
par.psi0=psi0;

% Utilities structure
u=@(c) CRRA(c,gamma);                              utility.u=u;
du=@(c) c.^(-gamma);                               utility.du=du;
duinv=@(du) du.^(-1/gamma);                        utility.duinv=duinv;
ddu=@(c) (-gamma)*c.^(-gamma-1);                   utility.ddu=ddu;
v=@(l) D*l.^(1+eta)./(1+eta);                      utility.v=v;
dv=@(l) D*l.^eta;                                  utility.dv=dv;
ddv=@(l) D*eta*l.^(eta-1);                         utility.ddv=ddv;
chi=@(tauL,tauLbar) (chi0/2)*(tauL-tauLbar).^2;    utility.chi=chi;
chi_tauL=@(tauL,tauLbar) chi0*(tauL-tauLbar);      utility.chi_tauL=chi_tauL;
chi_tauLbar=@(tauL,tauLbar) -chi0*(tauL-tauLbar);  utility.chi_tauLbar=chi_tauLbar;
psi=@(tauK,tauKbar) (psi0/2)*(tauK-tauKbar).^2;    utility.psi=psi;
psi_tauK=@(tauK,tauKbar) psi0*(tauK-tauKbar);      utility.psi_tauK=psi_tauK;
psi_tauKbar=@(tauK,tauKbar) -psi0*(tauK-tauKbar);  utility.psi_tauKbar=psi_tauKbar;

% Equation 15 David Stockman, get FC SS Rec. lag multiplier
lambda_FC = @(sigma) du(cFC)+sigma.*(du(cFC)+ddu(cFC).*cFC);
eq15 = @(sigma) -dv(lFC)+lambda_FC(sigma).*(1-alpha).*Z.*kFC.^(alpha).*lFC.^(-alpha)+sigma.*(-dv(lFC)-ddv(lFC).*lFC);
[sigmaFC,fvalFC]=fzero(eq15,0.0);

%% Generate Markov chain
[g_grd, Pg] = rouwen(gnum,log(g),rho,sigmasq);
g_grd=exp(g_grd);

%Generate shock sequences of TFP
chain = zeros(1,T);
chain(1)=1;
g_shocks=zeros(1,length(chain));
g_shocks_indeces=zeros(1,length(chain));
g_shocks(1)= g_grd(chain(1));
g_shocks_indeces(1)= chain(1);
for i=2:T
    this_step_distribution = Pg(chain(i-1),:);
    cumulative_distribution = cumsum(this_step_distribution);
    r = rand();
    chain(i) = find(cumulative_distribution>r,1);
    g_shocks_indeces(i)=chain(i);
    g_shocks(i)= g_grd(chain(i));
end

%% Countercheck parameters g process
log_g_stat=log(g_shocks)';
log_g_Y = log_g_stat(2:end);
log_g_X = log_g_stat(1:end-1)-log(g);

[beta_g, ~, res_g] = regress(log_g_Y, [ones(length(log_g_X),1) log_g_X]);
mu_g=exp(beta_g(1));
rho_g=beta_g(2);
sigma_inn_g = std(res_g);
%% Countercheck end

%% Projection parameters
order=3;
cross=true;
num_nodes=5;% #nodes in the grids
options=optimoptions(@fsolve,'Algorithm','Levenberg-Marquardt','MaxFunEvals',10000,'MaxIterations',10000,'Display','off','StepTolerance', 1.0000e-010);
options2=optimoptions(@fminunc,'Display','off','MaxFunEvals',10000,'MaxIterations',10000);

%% Solve FC stochastic (state contingent) dynamics
display(['*** FC solver! Chebyshev Order: ' num2str(order) ' Cross terms: ' num2str(cross) ' ***'])
display('Solving for FC stoch. dynamics...')
tic

Script_FCStochastic

% Simulate FC
k_FC(1)=kFC*1.00;
sigma_FC(1)=sigmaFC*1.0;
for t=1:T
    sigma_FC(t+1)=sigmatilde(k_FC(t),sigma_FC(t),g_shocks(t),phisigma_FC);
    c_FC(t)=ctilde(k_FC(t),sigma_FC(t),g_shocks(t),phic_FC);
    
    knextF=@(l) Z*k_FC(t)^alpha.*l.^(1-alpha)+(1-delta)*k_FC(t)-g_shocks(t)-c_FC(t);
    lambda_FCF = @(l) du(c_FC(t))+...
        -(sigma_FC(t+1)-sigma_FC(t)).*ddu(c_FC(t)).*knextF(l)+...
        +sigma_FC(t).*(du(c_FC(t))+ddu(c_FC(t)).*c_FC(t));
    eq15 = @(l) -dv(l)+lambda_FCF(l).*(1-alpha).*Z.*k_FC(t).^(alpha).*l.^(-alpha)+sigma_FC(t).*(-dv(l)-ddv(l).*l);
    [labor_FC(t),fval_FC(t)]=fzero(eq15,1.0);
    wage_FC(t)=(1-alpha)*Z*k_FC(t).^(alpha).*labor_FC(t).^(-alpha);
    tauL_FC(t)=1-(dv(labor_FC(t))/du(c_FC(t)))*(1/wage_FC(t));
    tauK_FC(t)=(g_shocks(t)-tauL_FC(t)*wage_FC(t)*labor_FC(t))/((alpha*Z*k_FC(t).^(alpha-1).*labor_FC(t).^(1-alpha)-delta)*k_FC(t));
    k_FC(t+1)=knextF(labor_FC(t));
end
Sequences{1}.g_shocks=g_shocks;Sequences{1}.c=c_FC;Sequences{1}.labor=labor_FC;Sequences{1}.k=k_FC;Sequences{1}.tauL=tauL_FC;Sequences{1}.tauK=tauK_FC;Sequences{1}.name='Stockman';

display(['FC mean - c: ' num2str(mean(c_FC(50:T))) ' K: ' num2str(mean(k_FC(50:T))) ' tauL: ' num2str(mean(tauL_FC(50:T))) ' tauK: ' num2str(mean(tauK_FC(50:T)))])
display(['FC std - c: ' num2str(std(c_FC(50:T))) ' K: ' num2str(std(k_FC(50:T))) ' tauL: ' num2str(std(tauL_FC(50:T))) ' tauK: ' num2str(std(tauK_FC(50:T)))])
time=toc;
display(['***********' num2str(time) 's ******************']);
display(' ');

%% Solve LTC stochastic Farhi dynamics with FOCs
display(['*** LTC FOCs Farhi solver! Chebyshev Order: ' num2str(order) ' Cross terms: ' num2str(cross) ' ***'])
display('Solving for LTC stoch. dynamics using FOCs...')
tic

Script_LTCStochastic_FOCs_Farhi

k_LTC_Fahri(1)=k_FC(1);
tauK_LTC_Fahri(1)=tauK_FC(1);
tauLbar_LTC_Fahri(1)=tauL_FC(1);
sigma_LTC_Fahri(1)=sigma_FC(1);
% Simulate LTC sequence
for t=1:T
    tauL_LTC_Fahri(t)=tauLtilde(k_LTC_Fahri(t),tauK_LTC_Fahri(t),g_shocks(t),sigma_LTC_Fahri(t),phitauL_LTC_Fahri);
    sigma_LTC_Fahri(t+1)=sigmatilde(k_LTC_Fahri(t),tauK_LTC_Fahri(t),g_shocks(t),sigma_LTC_Fahri(t),phisigma_LTC_Fahri);
    [labor_LTC_Fahri(t),c_LTC_Fahri(t),k_LTC_Fahri(t+1),labor_LTC_Fahri_tauL(t),c_LTC_Fahri_tauL(t),labor_LTC_Fahri_tauK(t),c_LTC_Fahri_tauK(t),labor_LTC_Fahri_k(t),c_LTC_Fahri_k(t)]=...
        getLaborCKnext(k_LTC_Fahri(t),tauK_LTC_Fahri(t),tauL_LTC_Fahri(t),Z,g_shocks(t),alpha,delta,gamma,eta,D);
    tauK_LTC_Fahri(t+1)=tauKtilde(k_LTC_Fahri(t),tauK_LTC_Fahri(t),g_shocks(t),sigma_LTC_Fahri(t),phitauK_LTC_Fahri);
    
    lambda_LTC_Fahri(t)=(du(c_LTC_Fahri(t))-(sigma_LTC_Fahri(t+1)-sigma_LTC_Fahri(t)).*ddu(c_LTC_Fahri(t)).*k_LTC_Fahri(t+1)+sigma_LTC_Fahri(t).*(du(c_LTC_Fahri(t))+ddu(c_LTC_Fahri(t)).*c_LTC_Fahri(t))+((-dv(labor_LTC_Fahri(t))+sigma_LTC_Fahri(t).*(-dv(labor_LTC_Fahri(t))-ddv(labor_LTC_Fahri(t)).*labor_LTC_Fahri(t)))).*labor_LTC_Fahri_tauL(t)./c_LTC_Fahri_tauL(t))./(1-(1-alpha).*Z.*k_LTC_Fahri(t).^(alpha).*labor_LTC_Fahri(t).^(-alpha).*labor_LTC_Fahri_tauL(t)./c_LTC_Fahri_tauL(t));
    nu_LTC_Fahri(t) = du(c_LTC_Fahri(t))-(sigma_LTC_Fahri(t+1)-sigma_LTC_Fahri(t)).*ddu(c_LTC_Fahri(t)).*k_LTC_Fahri(t+1)+sigma_LTC_Fahri(t).*(du(c_LTC_Fahri(t))+ddu(c_LTC_Fahri(t)).*c_LTC_Fahri(t))-lambda_LTC_Fahri(t);
    mu_LTC_Fahri(t)=(nu_LTC_Fahri(t).*c_LTC_Fahri_tauL(t))./labor_LTC_Fahri_tauL(t);
    
end
Sequences{2}.g_shocks=g_shocks;Sequences{2}.c=c_LTC_Fahri;Sequences{2}.labor=labor_LTC_Fahri;Sequences{2}.k=k_LTC_Fahri;Sequences{2}.tauL=tauL_LTC_Fahri;Sequences{2}.tauK=tauK_LTC_Fahri;Sequences{2}.name=['Fahri'];

display(['LTC Farhi mean - c: ' num2str(mean(c_LTC_Fahri(50:T))) ' K: ' num2str(mean(k_LTC_Fahri(50:T))) ' tauL: ' num2str(mean(tauL_LTC_Fahri(50:T))) ' tauK: ' num2str(mean(tauK_LTC_Fahri(50:T)))])
display(['LTC Farhi std - c: ' num2str(std(c_LTC_Fahri(50:T))) ' K: ' num2str(std(k_LTC_Fahri(50:T))) ' tauL: ' num2str(std(tauL_LTC_Fahri(50:T))) ' tauK: ' num2str(std(tauK_LTC_Fahri(50:T)))])

time=toc;
display(['***********' num2str(time) 's ******************']);
display(' ');


chi=@(tauL,tauLbar) (chi0/2)*(tauL-tauLbar).^2;    utility.chi=chi;
chi_tauL=@(tauL,tauLbar) chi0*(tauL-tauLbar);      utility.chi_tauL=chi_tauL;
chi_tauLbar=@(tauL,tauLbar) -chi0*(tauL-tauLbar);  utility.chi_tauLbar=chi_tauLbar;

save('./INIT/FCFARHI_temp')