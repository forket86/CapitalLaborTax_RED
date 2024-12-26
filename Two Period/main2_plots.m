%This file plots the results of the two period model


clear all
close all
clc

%load results workspace
load results2P

%set = 1 to save figures in /Figures folder, or = 0 not to save
save_figs = 1
figsSavePath = '../Figures/';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FC PLOTS

MM = 1; NN = 3;

figure('position',[10,10,800,225])

subplot(MM,NN,1)
plot(gammaKgrid,tauLFCc(:,1),'-','linewidth',2)
hold on
plot(gammaKgrid,tauLFCc(:,2),'--','linewidth',2)
grid on
xlabel('\gamma')
title('Labour tax: \tau^L(g)')
legend('low g','high g','Position',[0.255 0.710 0.082 0.117]);


subplot(MM,NN,2)
plot(gammaKgrid,tauKFCc(:,1),'-','linewidth',2)
hold on
plot(gammaKgrid,tauKFCc(:,2),'--','linewidth',2)
grid on
xlabel('\gamma')
title('Capital tax: \tau^K(g)')
legend('low g','high g')

subplot(MM,NN,3)
plot(gammaKgrid,tauKbarFCc,'-','linewidth',2)
hold on
plot(gammaKgrid,tauKFCc*Pg,'--','linewidth',2)
grid on
xlabel('\gamma')
legend('promise taubar','average tax','location','southeast')
title('Capital tax average and promise')

if save_figs
    savefig('../Figures/FigureB1a.fig')
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NC PLOTS

MM = 1; NN = 3;

figure('position',[10,10,800,225])

subplot(MM,NN,1)
plot(gammaKgrid,tauLNCc(:,1),'-','linewidth',2)
hold on
plot(gammaKgrid,tauLNCc(:,2),'--','linewidth',2)
grid on
xlabel('\gamma')
title('Labour tax: \tau^L(g)')
legend('low g','high g','location','southeast')

subplot(MM,NN,2)
plot(gammaKgrid,tauKNCc(:,1),'-','linewidth',2)
hold on
plot(gammaKgrid,tauKNCc(:,2),'--','linewidth',2)
grid on
xlabel('\gamma')
title('Capital tax: \tau^K(g)')
legend('low g','high g')

subplot(MM,NN,3)
plot(gammaKgrid,tauKbarNCc,'-','linewidth',2)
hold on
plot(gammaKgrid,tauKNCc*Pg,'--','linewidth',2)
grid on
xlabel('\gamma')
legend('promise taubar','average tax')
title('Capital tax average and promise')

if save_figs
    savefig('../Figures/FigureB1b.fig')
end




