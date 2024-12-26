T0=40;Tend=52;

time=-2:1:10;

hline_base = 3;
hline_comp = 1.5; 
fsize = 14;

%% Fig 10
figure
subplot(2,1,1)
hold on
plot(time,tauK_GEE(T0:Tend),'linewidth',hline_base)
plot(time,tauKbar_GEE(T0:Tend),'--','linewidth',hline_comp)
legend('$\tau^k$','$\overline{\tau}^k$', 'fontsize', fsize,'interpreter','latex')
ylabel('$\tau^k$', 'fontsize', fsize,'interpreter','latex')
grid on

%%
subplot(2,1,2)
hold on
plot(time,tauL_GEE(T0:Tend),'linewidth',hline_base)
plot(time,tauLbar_GEE(T0:Tend),'--','linewidth',hline_comp)
legend('$\tau^l$','$\overline{\tau}^l$', 'fontsize', fsize,'interpreter','latex')
ylabel('$\tau^l$', 'fontsize', fsize,'interpreter','latex')
grid on

savefig('./Figures/Figure10.fig')
close all

addpath('./LTC BB/');
rng(1)
infinityHorizon=10000;
T=infinityHorizon;
Tend=T;

SimulateGEE2

k_GEE=k(1:T);
c_GEE=c(1:T);
labor_GEE=labor(1:T);
tauL_GEE=tauL(1:T);
tauK_GEE=tauK(1:T);
tauLbar_GEE=tauLbar(1:T);
tauKbar_GEE=tauKbar(1:T);

U=u(c_GEE)-v(labor_GEE)-chi(tauL_GEE,tauLbar_GEE)-psi(tauK_GEE,tauKbar_GEE);
UtauL=u(c_GEE)-v(labor_GEE)-chi(tauL_GEE,tauLbar_GEE);
UtauK=u(c_GEE)-v(labor_GEE)-psi(tauK_GEE,tauKbar_GEE);
UHH=u(c_GEE)-v(labor_GEE);
V_GEE2=zeros(T-CUTEND_WELFARE,1);
VtauL_GEE=zeros(T-CUTEND_WELFARE,1);
VtauK_GEE=zeros(T-CUTEND_WELFARE,1);
VHH_GEE2=zeros(T-CUTEND_WELFARE,1);
for t=1:T-CUTEND_WELFARE
    gridGEE=0:length(U(t:t+CUTEND_WELFARE))-1;
    V_GEE2(t)=beta.^gridGEE*U(t:t+CUTEND_WELFARE)';
    VtauL_GEE(t)=beta.^gridGEE*UtauL(t:t+CUTEND_WELFARE)';
    VtauK_GEE(t)=beta.^gridGEE*UtauK(t:t+CUTEND_WELFARE)';
    VHH_GEE2(t)=beta.^gridGEE*UHH(t:t+CUTEND_WELFARE)';
end

GEEsamecost1_2=[
mean(tauK_GEE(T0:Tend));
mean(tauL_GEE(T0:Tend));
std(log(tauK_GEE(T0:Tend)+1));
std(log(tauL_GEE(T0:Tend)+1));
corr(log(tauK_GEE(T0+1:Tend)+1)',log(tauK_GEE(T0:Tend-1)+1)');
corr(log(tauL_GEE(T0+1:Tend)+1)',log(tauL_GEE(T0:Tend-1)+1)');
corr(log(tauK_GEE(T0:Tend)+1)',log(g_shocks(T0:Tend))');
corr(log(tauL_GEE(T0:Tend)+1)',log(g_shocks(T0:Tend))');];

Y=Z*k_GEE(1:Tend).^alpha.*labor_GEE(1:Tend).^(1-alpha);
I=k_GEE(2:Tend)-(1-delta)*k_GEE(1:Tend-1);
GEEsamecost2_2=[
mean(Y);
mean(c_GEE(T0:Tend));
mean(labor_GEE(T0:Tend));
mean(I);
std(Y);
std(log(c_GEE(T0:Tend)));
std(log(labor_GEE(T0:Tend)));
std(log(I));
corr(log(Y(T0+1:Tend))',log(Y(T0:Tend-1))');
corr(log(c_GEE(T0+1:Tend))',log(c_GEE(T0:Tend-1))');
corr(log(labor_GEE(T0+1:Tend))',log(labor_GEE(T0:Tend-1))');
corr(log(I(T0+1:end))',log(I(T0:end-1))');
corr(log(Y(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(c_GEE(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(labor_GEE(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(I(T0:end))',log(g_shocks(T0:Tend-1))');
];

GEEsamecost3_2=[
mean(tauKbar_GEE(T0:Tend));
mean(tauLbar_GEE(T0:Tend));
std(log(tauKbar_GEE(T0:Tend)+1));
std(log(tauLbar_GEE(T0:Tend)+1));
mean(tauL_GEE(T0:Tend)-tauLbar_GEE(T0:Tend));
mean(tauK_GEE(T0:Tend)-tauKbar_GEE(T0:Tend));
std(log(abs(tauL_GEE(T0:Tend)-tauLbar_GEE(T0:Tend))));
std(log(abs(tauK_GEE(T0:Tend)-tauKbar_GEE(T0:Tend))));];