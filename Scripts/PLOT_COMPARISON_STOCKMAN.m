%% Fig 1
%close all
T0=40;Tend=T0+12;

time=-2:1:10;

hline_base = 3;
hline_comp = 1.5; 
fsize = 14;

figure
subplot(3,1,1)
hold on
plot(time, g_shocks(T0:Tend), 'linewidth',hline_base)
ylabel('$g$', 'fontsize', fsize, 'interpreter','latex')

grid on
%%
subplot(3,1,2)
hold on
plot(time,tauK(T0:Tend),'linewidth',hline_base)
plot(time,tauK_FC(T0:Tend),'--','linewidth',hline_comp)
plot(time,tauK_LTC_Fahri(T0:Tend),'k-.','linewidth',hline_comp)
legend('Baseline','No CSC','Predet. $\tau^k$', 'fontsize', fsize,'interpreter','latex')
%legend('$\gamma_0^k=3.4,\  \gamma_0^l=140$','$\gamma_0^k=0.0,\  \gamma_0^l=0.0$','$\gamma_0^k=\infty,\  \gamma_0^l=0.0$', 'fontsize', fsize,'interpreter','latex')
ylabel('$\tau^k$', 'fontsize', fsize,'interpreter','latex')
grid on

%%
subplot(3,1,3)
hold on
plot(time,tauL(T0:Tend),'linewidth',hline_base)
plot(time,tauL_FC(T0:Tend),'--','linewidth',hline_comp)
plot(time,tauL_LTC_Fahri(T0:Tend),'k-.','linewidth',hline_comp)
ylabel('$\tau^l$', 'fontsize', fsize,'interpreter','latex')
grid on

savefig('./Figures/Figure1.fig')
close all

%% Fig 2
figure
subplot(2,2,1)
hold on
plot(time,labor(T0:Tend),'linewidth',hline_base)
plot(time,labor_FC(T0:Tend),'--','linewidth',hline_comp)
plot(time,labor_LTC_Fahri(T0:Tend),'k-.','linewidth',hline_comp)
ylabel('$l$', 'fontsize', fsize,'interpreter','latex')
grid on

%%
subplot(2,2,2)
hold on
plot(time,c(T0:Tend),'linewidth',hline_base)
plot(time,c_FC(T0:Tend),'--','linewidth',hline_comp)
plot(time,c_LTC_Fahri(T0:Tend),'k-.','linewidth',hline_comp)
ylabel('$c$', 'fontsize', fsize,'interpreter','latex')
legend('Baseline','No CSC','Predet. $\tau^k$', 'fontsize', fsize,'interpreter','latex')
%legend('$\gamma_0^k=3.4,\  \gamma_0^l=140$','$\gamma_0^k=0.0,\  \gamma_0^l=0.0$','$\gamma_0^k=\infty,\  \gamma_0^l=0.0$', 'fontsize', fsize,'interpreter','latex')
grid on

%%
subplot(2,2,3)
hold on
plot(time,k(T0:Tend),'linewidth',hline_base)
plot(time,k_FC(T0:Tend),'--','linewidth',hline_comp)
plot(time,k_LTC_Fahri(T0:Tend),'k-.','linewidth',hline_comp)
ylabel('$k$', 'fontsize', fsize,'interpreter','latex')
grid on

%%
subplot(2,2,4)
hold on
plot(time,Z*k(T0:Tend).^alpha.*labor(T0:Tend).^(1-alpha),'linewidth',hline_base)
plot(time,Z*k_FC(T0:Tend).^alpha.*labor_FC(T0:Tend).^(1-alpha),'--','linewidth',hline_comp)
plot(time,Z*k_LTC_Fahri(T0:Tend).^alpha.*labor_LTC_Fahri(T0:Tend).^(1-alpha),'k-.','linewidth',hline_comp)
ylabel('$y$', 'fontsize', fsize,'interpreter','latex')
grid on

savefig('./Figures/Figure2.fig')
close all

figure
%%
subplot(2,1,1)
hold on
plot(time,tauK(T0:Tend),'linewidth',hline_base)
plot(time,tauKbar(T0:Tend),'--','linewidth',hline_comp)
ylabel('$\tau^k$', 'fontsize', fsize,'interpreter','latex')
legend('$\tau^k$','$\overline{\tau}^k$', 'fontsize', fsize,'interpreter','latex')
grid on

%%
subplot(2,1,2)
hold on
plot(time,tauL(T0:Tend),'linewidth',hline_base)
plot(time,tauLbar(T0:Tend),'--','linewidth',hline_comp)
ylabel('$\tau^l$', 'fontsize', fsize,'interpreter','latex')
legend('$\tau^l$','$\overline{\tau}^l$', 'fontsize', fsize,'interpreter','latex')
grid on

savefig('./Figures/Figure3.fig')
close all

%% Generate Markov chain
addpath('./FC BB/');
rng(1)
infinityHorizon=10000;
T=infinityHorizon;

Simulate

tauLbar_temp = tauLbar(1:end-1);
tauKbar_temp = tauKbar(1:end-1);
U=u(c)-v(labor)-chi(tauL,tauLbar_temp)-psi(tauK,tauKbar_temp);
UtauL=u(c)-v(labor)-chi(tauL,tauLbar_temp);
UtauK=u(c)-v(labor)-psi(tauK,tauKbar_temp);
UHH=u(c)-v(labor);
V_FC=zeros(T-CUTEND_WELFARE,1);
VtauL_FC=zeros(T-CUTEND_WELFARE,1);
VtauK_FC=zeros(T-CUTEND_WELFARE,1);
VHH_FC=zeros(T-CUTEND_WELFARE,1);
for t=1:T-CUTEND_WELFARE
    gridC=0:length(U(t:t+CUTEND_WELFARE))-1;
    V_FC(t)=beta.^gridC*U(t:t+CUTEND_WELFARE)';
    VtauL_FC(t)=beta.^gridC*UtauL(t:t+CUTEND_WELFARE)';
    VtauK_FC(t)=beta.^gridC*UtauK(t:t+CUTEND_WELFARE)';
    VHH_FC(t)=beta.^gridC*UHH(t:t+CUTEND_WELFARE)';
end
T0=50;Tend=T;

FCsamecost1=[
mean(tauK(T0:Tend));
mean(tauL(T0:Tend));
std(log(tauK(T0:Tend)+1));
std(log(tauL(T0:Tend)+1));
corr(log(tauK(T0+1:Tend)+1)',log(tauK(T0:Tend-1)+1)');
corr(log(tauL(T0+1:Tend)+1)',log(tauL(T0:Tend-1)+1)');
corr(log(tauK(T0:Tend)+1)',log(g_shocks(T0:Tend))');
corr(log(tauL(T0:Tend)+1)',log(g_shocks(T0:Tend))');];

Y=Z*k(1:Tend).^alpha.*labor(1:Tend).^(1-alpha);
I=k(2:Tend)-(1-delta)*k(1:Tend-1);
FCsamecost2=[
mean(Y);
mean(c(T0:Tend));
mean(labor(T0:Tend));
mean(I);
std(Y);
std(log(c(T0:Tend)));
std(log(labor(T0:Tend)));
std(log(I));
corr(log(Y(T0+1:Tend))',log(Y(T0:Tend-1))');
corr(log(c(T0+1:Tend))',log(c(T0:Tend-1))');
corr(log(labor(T0+1:Tend))',log(labor(T0:Tend-1))');
corr(log(I(T0+1:end))',log(I(T0:end-1))');
corr(log(Y(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(c(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(labor(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(I(T0:end))',log(g_shocks(T0:Tend-1))');
];

FCsamecost3=[
mean(tauKbar(T0:Tend));
mean(tauLbar(T0:Tend));
std(log(tauKbar(T0:Tend)+1));
std(log(tauLbar(T0:Tend)+1));
mean(tauL(T0:Tend)-tauLbar(T0:Tend));
mean(tauK(T0:Tend)-tauKbar(T0:Tend));
std(log(abs(tauL(T0:Tend)-tauLbar(T0:Tend))));
std(log(abs(tauK(T0:Tend)-tauKbar(T0:Tend))));];

%%

rng(1)

[g_grd, Pg] = rouwen(gnum,log(g),rho,sigmasq);
g_grd=exp(g_grd);
T=10000;
Tend=T;

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

LOW_C=cFC*0.1;HIGH_C=cFC*1.9;
cFCunit=chebyTransformAle(LOW_C,HIGH_C,cFC,0);
LOW_K=kFC*(1-0.1);HIGH_K=kFC*(1+0.1);
LOW_TAUK=-0.4;HIGH_TAUK=0.4;
LOW_TAUL=tauLFC*(1-0.4); HIGH_TAUL=tauLFC*(1+0.4);
LOW_SIGMA=sigmaFC*(1-0.2);HIGH_SIGMA=sigmaFC*(1+0.2);
sigmaFCunit=chebyTransformAle(LOW_SIGMA,HIGH_SIGMA,sigmaFC,0);

% Grids structure
k_grd=linspace(LOW_K,HIGH_K,num_nodes)';           grids.k_grd=k_grd;
tauK_grd=linspace(LOW_TAUK,HIGH_TAUK,num_nodes)';  grids.tauK_grd=tauK_grd;
tauL_grd=linspace(LOW_TAUL,HIGH_TAUL,num_nodes)';  grids.tauL_grd=tauL_grd;
sigma_grd=linspace(LOW_SIGMA,HIGH_SIGMA,10);       grids.sigma_grd=sigma_grd;
grids.g_grd=g_grd; grids.z_grd=z_grd;
[grids.k_s,grids.sigma_s,grids.giter_s]=ndgrid(k_grd,sigma_grd,1:length(g_grd));

sigmatilde=@(k,sigma,g,phisigma) chebyAle(order,[k,sigma,g],phisigma,[LOW_K,LOW_SIGMA,g_grd(1)],[HIGH_K HIGH_SIGMA g_grd(end)],LOW_SIGMA,HIGH_SIGMA,cross);
ctilde    =@(k,sigma,g,phic)     chebyAle(order,[k,sigma,g],phic,    [LOW_K,LOW_SIGMA,g_grd(1)],[HIGH_K HIGH_SIGMA g_grd(end)],LOW_C,HIGH_C,cross);


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

CUTEND_WELFARE=500;
T0=50;Tend=T;
U=u(c_FC)-v(labor_FC);
V_STOCK=zeros(T-CUTEND_WELFARE,1);
for t=1:T-CUTEND_WELFARE
    gridB=0:length(U(t:t+CUTEND_WELFARE))-1;
    V_STOCK(t)=beta.^gridB*U(t:t+CUTEND_WELFARE)';
end

Stockman1=[
mean(tauK_FC(T0:Tend));
mean(tauL_FC(T0:Tend));
std(log(tauK_FC(T0:Tend)+1));
std(log(tauL_FC(T0:Tend)+1));
corr(log(tauK_FC(T0+1:Tend)+1)',log(tauK_FC(T0:Tend-1)+1)');
corr(log(tauL_FC(T0+1:Tend)+1)',log(tauL_FC(T0:Tend-1)+1)');
corr(log(tauK_FC(T0:Tend)+1)',log(g_shocks(T0:Tend))');
corr(log(tauL_FC(T0:Tend)+1)',log(g_shocks(T0:Tend))');];

Y=Z*k_FC(1:Tend).^alpha.*labor_FC(1:Tend).^(1-alpha);
I=k_FC(2:Tend)-(1-delta)*k_FC(1:Tend-1);
Stockman2=[
mean(Y);
mean(c_FC(T0:Tend));
mean(labor_FC(T0:Tend));
mean(I);
std(Y);
std(log(c_FC(T0:Tend)));
std(log(labor_FC(T0:Tend)));
std(log(I));
corr(log(Y(T0+1:Tend))',log(Y(T0:Tend-1))');
corr(log(c_FC(T0+1:Tend))',log(c_FC(T0:Tend-1))');
corr(log(labor_FC(T0+1:Tend))',log(labor_FC(T0:Tend-1))');
corr(log(I(T0+1:end))',log(I(T0:end-1))');
corr(log(Y(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(c_FC(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(labor_FC(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(I(T0:end))',log(g_shocks(T0:Tend-1))');
];

Stockman3=[NaN;NaN;NaN;NaN;NaN;NaN;NaN;NaN];
