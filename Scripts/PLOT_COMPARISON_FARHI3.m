% Simulation Column

load('./FC BB/INIT/chi0_39_psi0_39.mat')
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

load('./FC BB/INIT/chi0_39_psi0_0dot01.mat')

hline_base = 3;
hline_comp = 1.5; 
fsize = 14;
time=-2:1:10;

T0=40;Tend=52;
figure(4)
subplot(3,2,2)
hold on
plot(time, g_shocks(T0:Tend),'linewidth',hline_base)
ylabel('$g$','interpreter','latex')
grid on

%%
subplot(3,2,4)
hold on
plot(time,tauK100(T0:Tend),'linewidth',hline_base)
plot(time,tauK(T0:Tend),'--','linewidth',hline_comp)
legend('Baseline','$\gamma^k=0$', 'fontsize', fsize,'interpreter','latex')
ylabel('$\tau^K$','fontsize', fsize,'interpreter','latex')
grid on

%%
subplot(3,2,6)
hold on
plot(time,tauL100(T0:Tend),'linewidth',hline_base)
plot(time,tauL(T0:Tend),'--','linewidth',hline_comp)
ylabel('$\tau^L$','fontsize', fsize,'interpreter','latex')
grid on