%% Fig 5
T0=40;Tend=T0+12;

time=-2:1:10;

hline_base = 3;
hline_comp = 1.5; 
fsize = 14;

fig = figure(5)

left_color =  [0, 0.4470, 0.7410];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

subplot(2,2,1)
hold on
plot(time, g_shocks(T0:Tend), 'linewidth',hline_base)
ylabel('$g$', 'fontsize', fsize, 'interpreter','latex')

grid on
%
subplot(2,2,2)
hold on
yyaxis left
plot(time,tauK(T0:Tend),'linewidth',hline_base)
ylabel('$\tau^k$', 'fontsize', fsize,'interpreter','latex')
yyaxis right
haxes =plot(time,tauK_LTC_Fahri_Debt(T0:Tend),'k-.','linewidth',hline_comp)
legend('Baseline','Predet. $\tau^k$', 'fontsize', fsize,'interpreter','latex')
ylabel('$\tau^k$', 'fontsize', fsize,'interpreter','latex')
grid on

%
subplot(2,2,3)
hold on
plot(time,tauL(T0:Tend),'linewidth',hline_base)
plot(time,tauL_LTC_Fahri_Debt(T0:Tend),'k-.','linewidth',hline_comp)
ylabel('$\tau^l$', 'fontsize', fsize,'interpreter','latex')
grid on

%
subplot(2,2,4)
hold on
plot(time,b(T0:Tend),'linewidth',hline_base)
plot(time,b_LTC_Fahri_Debt(T0:Tend),'k-.','linewidth',hline_comp)
ylabel('$b$', 'fontsize', fsize,'interpreter','latex')
grid on