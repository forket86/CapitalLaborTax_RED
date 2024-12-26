clear all;
clc;
close all;

addpath('../Utils/');

clear all;
clc;

load('./Measurability Neural Integrand/Good/chi0_39_psi0_39a2.mat')

gaptauK=tauK-tauKbar(1:end-1);
gaptauL=tauL-tauLbar(1:end-1);

epsilon=0.01;
mean(gaptauK(find(gaptauK> epsilon)))
mean(gaptauK(find(gaptauK< -epsilon)))

mean(gaptauL(find(gaptauL> epsilon)))   
epsilon=0.023;
mean(gaptauL(find(gaptauL< -epsilon))) 

%% Fig 1
%close all
T0=40;Tend=52;

time=-2:1:10;

hline_base = 3;
hline_comp = 1.5; 
fsize = 14;

figure
subplot(3,2,1)
hold on
plot(time, g_shocks(T0:Tend), 'linewidth',hline_base)
xlim([-1 10])
ylabel('$g$', 'fontsize', fsize, 'interpreter','latex')

grid on
%
subplot(3,2,3)
hold on
plot(time,tauK(T0:Tend),'linewidth',hline_base)
plot(time,tauKbar(T0:Tend),'linewidth',hline_base)
ylim([-0.05 0.05])
xlim([-1 10])
ylabel('$\tau^k$', 'fontsize', fsize,'interpreter','latex')
grid on

%
subplot(3,2,5)
hold on
plot(time,tauL(T0:Tend),'linewidth',hline_base)
plot(time,tauLbar(T0:Tend),'linewidth',hline_base)
ylim([.28 0.35])
xlim([-1 10])
ylabel('$\tau^l$', 'fontsize', fsize,'interpreter','latex')
grid on

%%
T0=40+210;Tend=52+210;
time=-2:1:10;
subplot(3,2,2)
hold on
plot(time, g_shocks(T0:Tend), 'linewidth',hline_base)
xlim([-1 10])
ylabel('$g$', 'fontsize', fsize, 'interpreter','latex')
grid on

subplot(3,2,4)
hold on
plot(time,tauK(T0:Tend),'linewidth',hline_base)
plot(time,tauKbar(T0:Tend),'linewidth',hline_base)
ylim([-0.05 0.05])
xlim([-1 10])
ylabel('$\tau^k$', 'fontsize', fsize,'interpreter','latex')
grid on

%
subplot(3,2,6)
hold on
plot(time,tauL(T0:Tend),'linewidth',hline_base)
plot(time,tauLbar(T0:Tend),'linewidth',hline_base)
ylim([.28 0.35])
xlim([-1 10])
ylabel('$\tau^l$', 'fontsize', fsize,'interpreter','latex')
grid on
