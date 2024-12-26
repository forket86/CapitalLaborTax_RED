%% Table
warning off
load('./FC/INIT/FCFARHI.mat')
warning on

rng(1)

explore=0.1;
LOW_C=cFC*0.01;HIGH_C=cFC*1.5;
cFCunit=chebyTransformAle(LOW_C,HIGH_C,cFC,0);
LOW_K=kFC*(1-explore);HIGH_K=kFC*(1+explore);
LOW_TAUK=-0.15;HIGH_TAUK=0.15;
tauKFCunit=chebyTransformAle(LOW_TAUK,HIGH_TAUK,tauKFC,0);
LOW_TAUL=tauLFC*(1-explore); HIGH_TAUL=tauLFC*(1+explore);
tauLFCunit=chebyTransformAle(LOW_TAUL,HIGH_TAUL,tauLFC,0);
LOW_SIGMA=sigmaFC*(1-0.15);HIGH_SIGMA=sigmaFC*(1+0.1);
sigmaFCunit=chebyTransformAle(LOW_SIGMA,HIGH_SIGMA,sigmaFC,0);

% Grids structure
k_grd=linspace(LOW_K,HIGH_K,num_nodes)';           grids.k_grd=k_grd;
tauK_grd=linspace(LOW_TAUK,HIGH_TAUK,num_nodes)';  grids.tauK_grd=tauK_grd;
tauL_grd=linspace(LOW_TAUL,HIGH_TAUL,num_nodes)';  grids.tauL_grd=tauL_grd;
sigma_grd=linspace(LOW_SIGMA,HIGH_SIGMA,10);       grids.sigma_grd=sigma_grd;
grids.g_grd=g_grd; grids.z_grd=z_grd;
[grids.k_s,grids.tauK_s,grids.giter_s,grids.sigma_s]=ndgrid(k_grd,tauK_grd,1:length(g_grd),sigma_grd);

sigmatilde   =@(k,tauK,g,sigma,phisigma)    chebyAle(order,[k,tauK,g,sigma],phisigma,   [LOW_K,LOW_TAUK,g_grd(1),LOW_SIGMA],[HIGH_K,HIGH_TAUK,g_grd(end),HIGH_SIGMA],LOW_SIGMA,HIGH_SIGMA,cross);
tauKtilde    =@(k,tauK,g,sigma,phitauK)     chebyAle(order,[k,tauK,g,sigma],phitauK,    [LOW_K,LOW_TAUK,g_grd(1),LOW_SIGMA],[HIGH_K,HIGH_TAUK,g_grd(end),HIGH_SIGMA],LOW_TAUK,HIGH_TAUK,cross);
tauLtilde    =@(k,tauK,g,sigma,phitauL)     chebyAle(order,[k,tauK,g,sigma],phitauL,    [LOW_K,LOW_TAUK,g_grd(1),LOW_SIGMA],[HIGH_K,HIGH_TAUK,g_grd(end),HIGH_SIGMA],LOW_TAUL,HIGH_TAUL,cross);

T0=50;Tend=10000;

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
for i=2:Tend
    this_step_distribution = Pg(chain(i-1),:);
    cumulative_distribution = cumsum(this_step_distribution);
    r = rand();
    chain(i) = find(cumulative_distribution>r,1);
    g_shocks_indeces(i)=chain(i);
    g_shocks(i)= g_grd(chain(i));
end

k_LTC_Fahri(1)=k_FC(1);
tauK_LTC_Fahri(1)=tauK_FC(1);
tauLbar_LTC_Fahri(1)=tauL_FC(1);
sigma_LTC_Fahri(1)=sigma_FC(1);
% Simulate LTC sequence
for t=1:Tend
    tauL_LTC_Fahri(t)=tauLtilde(k_LTC_Fahri(t),tauK_LTC_Fahri(t),g_shocks(t),sigma_LTC_Fahri(t),phitauL_LTC_Fahri);
    sigma_LTC_Fahri(t+1)=sigmatilde(k_LTC_Fahri(t),tauK_LTC_Fahri(t),g_shocks(t),sigma_LTC_Fahri(t),phisigma_LTC_Fahri);
    [labor_LTC_Fahri(t),c_LTC_Fahri(t),k_LTC_Fahri(t+1),labor_LTC_Fahri_tauL(t),c_LTC_Fahri_tauL(t),labor_LTC_Fahri_tauK(t),c_LTC_Fahri_tauK(t),labor_LTC_Fahri_k(t),c_LTC_Fahri_k(t)]=...
        getLaborCKnext(k_LTC_Fahri(t),tauK_LTC_Fahri(t),tauL_LTC_Fahri(t),Z,g_shocks(t),alpha,delta,gamma,eta,D);
    tauK_LTC_Fahri(t+1)=tauKtilde(k_LTC_Fahri(t),tauK_LTC_Fahri(t),g_shocks(t),sigma_LTC_Fahri(t),phitauK_LTC_Fahri);
    
    lambda_LTC_Fahri(t)=(du(c_LTC_Fahri(t))-(sigma_LTC_Fahri(t+1)-sigma_LTC_Fahri(t)).*ddu(c_LTC_Fahri(t)).*k_LTC_Fahri(t+1)+sigma_LTC_Fahri(t).*(du(c_LTC_Fahri(t))+ddu(c_LTC_Fahri(t)).*c_LTC_Fahri(t))+((-dv(labor_LTC_Fahri(t))+sigma_LTC_Fahri(t).*(-dv(labor_LTC_Fahri(t))-ddv(labor_LTC_Fahri(t)).*labor_LTC_Fahri(t)))).*labor_LTC_Fahri_tauL(t)./c_LTC_Fahri_tauL(t))./(1-(1-alpha).*Z.*k_LTC_Fahri(t).^(alpha).*labor_LTC_Fahri(t).^(-alpha).*labor_LTC_Fahri_tauL(t)./c_LTC_Fahri_tauL(t));
    nu_LTC_Fahri(t) = du(c_LTC_Fahri(t))-(sigma_LTC_Fahri(t+1)-sigma_LTC_Fahri(t)).*ddu(c_LTC_Fahri(t)).*k_LTC_Fahri(t+1)+sigma_LTC_Fahri(t).*(du(c_LTC_Fahri(t))+ddu(c_LTC_Fahri(t)).*c_LTC_Fahri(t))-lambda_LTC_Fahri(t);
    mu_LTC_Fahri(t)=(nu_LTC_Fahri(t).*c_LTC_Fahri_tauL(t))./labor_LTC_Fahri_tauL(t);
    
end


CUTEND_WELFARE=500;
U=u(c_LTC_Fahri)-v(labor_LTC_Fahri);
V_FAHRI=zeros(Tend-CUTEND_WELFARE,1);
for t=1:Tend-CUTEND_WELFARE
    grid=0:length(U(t:t+CUTEND_WELFARE))-1;
    V_FAHRI(t)=beta.^grid*U(t:t+CUTEND_WELFARE)';
end
Fahri1=[
mean(tauK_LTC_Fahri(T0:Tend));
mean(tauL_LTC_Fahri(T0:Tend));
std(log(tauK_LTC_Fahri(T0:Tend)+1));
std(log(tauL_LTC_Fahri(T0:Tend)+1));
corr(log(tauK_LTC_Fahri(T0+1:Tend)+1)',log(tauK_LTC_Fahri(T0:Tend-1)+1)');
corr(log(tauL_LTC_Fahri(T0+1:Tend)+1)',log(tauL_LTC_Fahri(T0:Tend-1)+1)');
corr(log(tauK_LTC_Fahri(T0:Tend)+1)',log(g_shocks(T0:Tend))');
corr(log(tauL_LTC_Fahri(T0:Tend)+1)',log(g_shocks(T0:Tend))')];

Y=Z*k_LTC_Fahri(1:Tend).^alpha.*labor_LTC_Fahri(1:Tend).^(1-alpha);
I=k_LTC_Fahri(2:Tend)-(1-delta)*k_LTC_Fahri(1:Tend-1);
Fahri2=[
mean(Y);
mean(c_LTC_Fahri(T0:Tend));
mean(labor_LTC_Fahri(T0:Tend));
mean(I);
std(Y);
std(log(c_LTC_Fahri(T0:Tend)));
std(log(labor_LTC_Fahri(T0:Tend)));
std(log(I));
corr(log(Y(T0+1:Tend))',log(Y(T0:Tend-1))');
corr(log(c_LTC_Fahri(T0+1:Tend))',log(c_LTC_Fahri(T0:Tend-1))');
corr(log(labor_LTC_Fahri(T0+1:Tend))',log(labor_LTC_Fahri(T0:Tend-1))');
corr(log(I(T0+1:end))',log(I(T0:end-1))');
corr(log(Y(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(c_LTC_Fahri(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(labor_LTC_Fahri(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(I(T0:end))',log(g_shocks(T0:Tend-1))');
];

Fahri3=[NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];