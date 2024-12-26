% Simulation Column
warning off
load('./FC BB/INIT/chi0_39_psi0_39.mat')
warning on
tauK100=tauK;
tauL100=tauL;
k100=k;
tauKbar100=tauKbar;
tauLbar100=tauLbar;
labor100=labor;
c100=c;

chi01=chi0;
psi01=psi0;

tauKF=tauK_LTC_Fahri;
tauLF=tauL_LTC_Fahri;
kF=k_LTC_Fahri;
laborF=labor_LTC_Fahri;
cF=c_LTC_Fahri;
warning off
load('./FC BB/INIT/chi0_0dot01_psi0_39.mat')
warning on
hline_base = 3;
hline_comp = 1.5; 
fsize = 14;
time=-2:1:10;

T0=40;Tend=52;
figure(4)
subplot(3,2,1)
hold on
plot(time, g_shocks(T0:Tend),'linewidth',hline_base)
ylabel('$g$','interpreter','latex')
grid on

%%
subplot(3,2,3)
hold on
plot(time,tauK100(T0:Tend),'linewidth',hline_base)
plot(time,tauK(T0:Tend),'--','linewidth',hline_comp)
legend('Baseline','$\gamma^l=0$', 'fontsize', fsize,'interpreter','latex')
ylabel('$\tau^K$','fontsize', fsize,'interpreter','latex')
%legend(['FC \chi_0=' num2str(chi01) '\psi_0=' num2str(psi01)],['FC \chi_0=' num2str(chi0) '\psi_0=' num2str(psi0)], 'Fahri')
grid on

%%
subplot(3,2,5)
hold on
plot(time,tauL100(T0:Tend),'linewidth',hline_base)
plot(time,tauL(T0:Tend),'--','linewidth',hline_comp)
%plot(time,tauLF(T0:Tend),'k-.','linewidth',hline_comp)
ylabel('$\tau^L$','fontsize', fsize,'interpreter','latex')
grid on