T0=40;Tend=52;

time=-2:1:10;

hline_base = 3;
hline_comp = 1.5; 
fsize = 14;

%% Fig 7
figure
subplot(3,1,1)
hold on
plot(time,g_shocks(T0:Tend),'linewidth',hline_base)
ylabel('$g$', 'fontsize', fsize,'interpreter','latex')
grid on

%%
subplot(3,1,2)
hold on
plot(time,tauK_GEE(T0:Tend),'linewidth',hline_base)
%plot(time,tauK_FCM(T0:Tend),'--','linewidth',hline_comp)
ylabel('$\tau^k$', 'fontsize', fsize,'interpreter','latex')
grid on

%%
subplot(3,1,3)
hold on
plot(time,tauL_GEE(T0:Tend),'linewidth',hline_base)
%plot(time,tauL_FCM(T0:Tend),'--','linewidth',hline_comp)
ylabel('$\tau^l$', 'fontsize', fsize,'interpreter','latex')
grid on

savefig('./Figures/Figure7.fig')
close all


%% Fig 8
figure

%%
subplot(2,2,1)
hold on
plot(time,labor_GEE(T0:Tend),'linewidth',hline_base)
%plot(time,labor_FCM(T0:Tend),'--','linewidth',hline_comp)
ylabel('$l$', 'fontsize', fsize,'interpreter','latex')
grid on

%%
subplot(2,2,2)
hold on
plot(time,c_GEE(T0:Tend),'linewidth',hline_base)
%plot(time,c_FCM(T0:Tend),'--','linewidth',hline_comp)
ylabel('$c$', 'fontsize', fsize,'interpreter','latex')
grid on

subplot(2,2,3)
hold on
plot(time,k_GEE(T0:Tend),'linewidth',hline_base)
%plot(time,k_FCM(T0:Tend),'--','linewidth',hline_comp)
grid on
ylabel('$k$', 'fontsize', fsize,'interpreter','latex')
%legend(['FC \chi_0=' num2str(chi0) '\psi_0=' num2str(psi0)],['GEE \chi_0=' num2str(chi0) '\psi_0=' num2str(psi0)])
%legend(['FC \chi_0=' num2str(chi0) '\psi_0=' num2str(psi0)],['GEE \chi_0=' num2str(chi0) '\psi_0=' num2str(psi0)],'NC')
%legend(['GEE \chi_0=' num2str(chi0) '\psi_0=' num2str(psi0)],'NC')

%%
subplot(2,2,4)
hold on
plot(time,Z*k_GEE(T0:Tend).^alpha.*labor_GEE(T0:Tend).^(1-alpha),'linewidth',hline_base)
%plot(time,Z*k_FCM(T0:Tend).^alpha.*labor_FCM(T0:Tend).^(1-alpha),'--','linewidth',hline_comp)
ylabel('$y$', 'fontsize', fsize,'interpreter','latex')
grid on

savefig('./Figures/Figure8.fig')
close all

%% Fig 9
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

savefig('./Figures/Figure9.fig')
close all


addpath('./LTC BB/');
rng(1)
infinityHorizon=10000;
T=infinityHorizon;
Tend=T;

SimulateGEE

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
V_GEE=zeros(T-CUTEND_WELFARE,1);
VtauL_GEE=zeros(T-CUTEND_WELFARE,1);
VtauK_GEE=zeros(T-CUTEND_WELFARE,1);
VHH_GEE=zeros(T-CUTEND_WELFARE,1);
for t=1:T-CUTEND_WELFARE
    gridGEE=0:length(U(t:t+CUTEND_WELFARE))-1;
    V_GEE(t)=beta.^gridGEE*U(t:t+CUTEND_WELFARE)';
    VtauL_GEE(t)=beta.^gridGEE*UtauL(t:t+CUTEND_WELFARE)';
    VtauK_GEE(t)=beta.^gridGEE*UtauK(t:t+CUTEND_WELFARE)';
    VHH_GEE(t)=beta.^gridGEE*UHH(t:t+CUTEND_WELFARE)';
end

GEEsamecost1=[
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
GEEsamecost2=[
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

GEEsamecost3=[
mean(tauKbar_GEE(T0:Tend));
mean(tauLbar_GEE(T0:Tend));
std(log(tauKbar_GEE(T0:Tend)+1));
std(log(tauLbar_GEE(T0:Tend)+1));
mean(tauL_GEE(T0:Tend)-tauLbar_GEE(T0:Tend));
mean(tauK_GEE(T0:Tend)-tauKbar_GEE(T0:Tend));
std(log(abs(tauL_GEE(T0:Tend)-tauLbar_GEE(T0:Tend))));
std(log(abs(tauK_GEE(T0:Tend)-tauKbar_GEE(T0:Tend))));];