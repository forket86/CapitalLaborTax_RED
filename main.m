%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file can be run to replicate all figures and tables in the paper.
% See Readme file for details.


clear all;
clc;
close all;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUANTITATIVE MODEL PLOTS:


addpath('./Utils/');
addpath('./Functions/');
addpath('./Scripts/');

%% Figures 4 and 5
warning off
load('./FC Debt - Farhi/INIT/035_Farhi_Extreme')
warning on
tauK_LTC_Fahri_Debt=tauK;
tauL_LTC_Fahri_Debt=tauL;
labor_LTC_Fahri_Debt=labor;
c_LTC_Fahri_Debt=c;
k_LTC_Fahri_Debt=k;
b_LTC_Fahri_Debt=b;
MomentsFarhi
warning off
load('./FC Debt/INIT/openbound35_2')
warning on
% Figure 5
PLOT_COMPARISON_DEBT
savefig('./Figures/Figure5.fig')
close all;

MomentsFCdebt

PLOT_COMPARISON_FARHI

% Figure 4
PLOT_COMPARISON_FARHI2
PLOT_COMPARISON_FARHI3
savefig('./Figures/Figure4.fig')
close all;

%% Fig 1, 2, and 3 comparison with stockman, same costs on both
%clear all;
clc;

warning off
load('./FC BB/INIT/chi0_39_psi0_39.mat')
warning on
PLOT_COMPARISON_STOCKMAN

%% Fig 7 and 8 GEE con same costs on both, compared with our FC
clc;
warning off
load('./FC BB/INIT/chi0_39_psi0_39.mat')
warning on
k_FCM=k;
c_FCM=c;
labor_FCM=labor;
tauL_FCM=tauL;
tauK_FCM=tauK;
tauLbar_FCM=tauLbar;
tauKbar_FCM=tauKbar;

warning off
load('./FC BB/INIT/CurveGEE/Both/Curve39')
warning on

k_GEE=k(1:T);
c_GEE=c(1:T);
labor_GEE=labor(1:T);
tauL_GEE=tauL(1:T);
tauK_GEE=tauK(1:T);
tauLbar_GEE=tauLbar(1:T);
tauKbar_GEE=tauKbar(1:T);
close all
PLOT_COMPARISON_GEE

tauKStore=tauK;
tauLStore=tauL;
sigmanextStore=sigma;
tauKbarnextStore=tauKbar;
tauLbarnextStore=tauLbar;

warning off
load('./FC BB/INIT/CurveGEE/Both/Curve2')
warning on

k_GEE=k(1:T);
c_GEE=c(1:T);
labor_GEE=labor(1:T);
tauL_GEE=tauL(1:T);
tauK_GEE=tauK(1:T);
tauLbar_GEE=tauLbar(1:T);
tauKbar_GEE=tauKbar(1:T);
close all


PLOT_COMPARISON_GEE2

rows1={'tauK';'tauL';'std(log(tauK+1))';'std(log(tauL+1))';'autocorr(log(tauK+1))';'autocorr(log(tauL+1))'};
rows2={'y';'c';'labor';'investment';'std(log(y))';'std(log(c))';'std(log(labor))';'std(log(investment))';'autocorr(log(y))';'autocorr(log(c))';'autocorr(log(labor))';'autocorr(log(investment))'};

labels={'Moment','BB CSC','BB no CSC','BB pred. tauK', 'Debt CSC', 'Debt pred. tauK','Data'}
%Table 1
data=[0.355
    0.226
    0.022
    0.015
    0.868
    0.876];

t1=table(rows1,round(FCsamecost1(1:6),3),round(Stockman1(1:6),3),round(Fahri1(1:6),3),round(FCdebtsamecost1(1:6),3),round(Fahrisamecost1(1:6),3),round(data,3),'VariableNames',labels)
table2latex(t1, './Tables/Table2.tex')

labels={'Moment','LTC (gamma=39)','LTC (gamma=2)','FC (gamma=39)','Data'}
t5=table(rows1,round(GEEsamecost1(1:6),3),round(GEEsamecost1_2(1:6),3),round(FCsamecost1(1:6),3),round(data,3),'VariableNames',labels)
table2latex(t5, './Tables/Table5.tex')

labels={'Moment','FC','FC no CSC','FC pred. tauK', 'LTC'}
t2=table(rows2,round(FCsamecost2(1:12),3),round(Stockman2(1:12),3),round(Fahri2(1:12),3),round(GEEsamecost2(1:12),3),'VariableNames',labels)
table2latex(t2, './Tables/TableC1.tex')



%% Table3_gAndTFPShocks
warning off
load('./FC BB/INIT/TFP39.mat')
warning on
T0=50;Tend=T;
FC39=[
std(log(tauK(T0:Tend)+1));
std(log(tauL(T0:Tend)+1));
corr(log(tauK(T0+1:Tend)+1)',log(tauK(T0:Tend-1)+1)');
corr(log(tauL(T0+1:Tend)+1)',log(tauL(T0:Tend-1)+1)');];

warning off
load('./FC BB/INIT/TFP53.mat')
warning on
T0=50;Tend=T;
FC53=[
std(log(tauK(T0:Tend)+1));
std(log(tauL(T0:Tend)+1));
corr(log(tauK(T0+1:Tend)+1)',log(tauK(T0:Tend-1)+1)');
corr(log(tauL(T0+1:Tend)+1)',log(tauL(T0:Tend-1)+1)');];

warning off
load('./FC BB/INIT/TFP0dot01_39_Fahri.mat')
warning on
T0=50;Tend=T;
FC0dot01_39=[
std(log(tauK(T0:Tend)+1));
std(log(tauL(T0:Tend)+1));
corr(log(tauK(T0+1:Tend)+1)',log(tauK(T0:Tend-1)+1)');
corr(log(tauL(T0+1:Tend)+1)',log(tauL(T0:Tend-1)+1)');];

rows1={'std(log(tauK+1))';'std(log(tauL+1))';'autocorr(log(tauK+1))';'autocorr(log(tauL+1))'};
labels={'Moment','FC TFP 53','FC TFP 39','FC pred. tauK'}
t=table(rows1,round(FC53,3),round(FC39,3),round(FC0dot01_39,3),'VariableNames',labels)
table2latex(t, './Tables/Table3.tex')

%% Table4_AsymmetricCost
warning off
load('./FC BB/INIT/chi0_39_psi0_39a2.mat')
warning on
T0=50;Tend=T;
FC39a=[
mean(tauK(T0:Tend));
mean(tauL(T0:Tend));
std(log(tauK(T0:Tend)+1));
std(log(tauL(T0:Tend)+1));
corr(log(tauK(T0+1:Tend)+1)',log(tauK(T0:Tend-1)+1)');
corr(log(tauL(T0+1:Tend)+1)',log(tauL(T0:Tend-1)+1)');];

rows1={'E(tauK)';'E(tauL)';'std(log(tauK+1))';'std(log(tauL+1))';'autocorr(log(tauK+1))';'autocorr(log(tauL+1))'};
labels={'Moment','Baseline (k=1)','Baseline (k=2)'}
t=table(rows1,round(FCsamecost1(1:6),3),round(FC39a,3),'VariableNames',labels)
table2latex(t, './Tables/Table4.tex')

%% TableC2_LowPersistenceOfg
warning off
load('./FC BB/INIT/LowPers39_39.mat')
warning on
T0=50;Tend=T;
Baseline39=[
std(log(tauK(T0:Tend)+1));
std(log(tauL(T0:Tend)+1));
corr(log(tauK(T0+1:Tend)+1)',log(tauK(T0:Tend-1)+1)');
corr(log(tauL(T0+1:Tend)+1)',log(tauL(T0:Tend-1)+1)');];

T0=50;Tend=T;
LTC39=[
std(log(tauK_LTC_Fahri(T0:Tend)+1));
std(log(tauL_LTC_Fahri(T0:Tend)+1));
corr(log(tauK_LTC_Fahri(T0+1:Tend)+1)',log(tauK_LTC_Fahri(T0:Tend-1)+1)');
corr(log(tauL_LTC_Fahri(T0+1:Tend)+1)',log(tauL_LTC_Fahri(T0:Tend-1)+1)');];

T0=50;Tend=T;
FC39=[
std(log(tauK_FC(T0:Tend)+1));
std(log(tauL_FC(T0:Tend)+1));
corr(log(tauK_FC(T0+1:Tend)+1)',log(tauK_FC(T0:Tend-1)+1)');
corr(log(tauL_FC(T0+1:Tend)+1)',log(tauL_FC(T0:Tend-1)+1)');];

rows1={'std(log(tauK+1))';'std(log(tauL+1))';'autocorr(log(tauK+1))';'autocorr(log(tauL+1))'};
labels={'Moment','BB CSC','BB no CSC','BB pred. \tau^k'}
t=table(rows1,round(Baseline39,3),round(FC39,3),round(LTC39,3),'VariableNames',labels)
table2latex(t, './Tables/TableC2.tex')

% Fig C2
T0=251;Tend=T0+12;

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

savefig('./Figures/FigureC2.fig')
close all

% Fig C3
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

savefig('./Figures/FigureC3.fig')
close all

%% Figure C1
warning off
% Set this to 1 to recalculate moments from scratch
RECALCULATE_SIMULATION=0;
infinityHorizon=5000;

if RECALCULATE_SIMULATION
    
    addpath('./FC BB/');
    
    cutInit=50;
    costgrid=[0.1:0.1:0.9 5 20 30 39 50];
    for costIndex=1:length(costgrid)
        
        if costgrid(costIndex)<1
            load(['./FC BB/INIT/Curve/Both/Curve0dot' num2str(costIndex)])
        else
            load(['./FC BB/INIT/Curve/Both/Curve' num2str(costgrid(costIndex))])
        end
        
        T=infinityHorizon;
        rng(1)

        Simulate
        
        ETauL(costIndex)=mean(tauL(cutInit:T));
        ETauK(costIndex)=mean(tauK(cutInit:T));
        StdTauL(costIndex)=std(log(tauL(cutInit:T)+1));
        StdTauK(costIndex)=std(log(tauK(cutInit:T)+1));
        
        if costgrid(costIndex)<1
            load(['./FC BB/INIT/Curve/Lowchi/CurveLOWCHI0dot' num2str(costIndex)])
        else
            load(['./FC BB/INIT/Curve/Lowchi/CurveLOWCHI' num2str(costgrid(costIndex))])
        end
        
        T=infinityHorizon;
        rng(1)

        Simulate
        
        ETauL2(costIndex)=mean(tauL(cutInit:T));
        ETauK2(costIndex)=mean(tauK(cutInit:T));
        StdTauL2(costIndex)=std(log(tauL(cutInit:T)+1));
        StdTauK2(costIndex)=std(log(tauK(cutInit:T)+1));
    end
else
    load('./FC BB/Simulations')
end
warning on

% figure
figure
hline_base = 3;
hline_data = 1;
hline_comp = 2;
fsize = 14;


subplot(2,1,1)
hold on
plot([0.01 costgrid],[std(log(tauK_FC(cutInit:end)+1)) StdTauK],'linewidth',hline_base)
plot([0.01 costgrid],[std(log(tauK_FC(cutInit:end)+1)) StdTauK2],'--','linewidth',hline_base)
plot([0.01 costgrid],0.0212*ones(1,length(costgrid)+1),'k-.','linewidth',hline_data)
title('$\sigma(\tau^k)$', 'fontsize', fsize,'interpreter','latex')
xlabel('$\gamma$', 'fontsize', fsize,'interpreter','latex')
grid on
xlim([0,50])
ylim([0,0.07])

subplot(2,1,2)
hold on
plot([0.01 costgrid],[std(log(tauL_FC(cutInit:end)+1)) StdTauL],'linewidth',hline_base)
plot([0.01 costgrid],[std(log(tauL_FC(cutInit:end)+1)) StdTauL2],'--','linewidth',hline_base)
plot([0.01 costgrid],0.0159*ones(1,length(costgrid)+1),'k-.','linewidth',hline_data)
title('$\sigma(\tau^l)$', 'fontsize', fsize,'interpreter','latex')
xlabel('$\gamma$', 'fontsize', fsize,'interpreter','latex')
grid on
xlim([0,50])

savefig('./Figures/FigureC1.fig')
close all;
clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA PLOTS:

% Figure D1: Tax rates and personal income tax counterfactuals in the data
run('./Data/main1_dataprep.m')
close all;
clc;

% Figures 6 and D2: Local Projection results
run('./Data/main3_plotLP.m')
close all;
clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TWO PERIOD MODEL PLOTS:

% Figure B1: Solution of two period model in FC and LTC cases
run('./Two Period/main2_plots.m')
close all;
clc;

display(['All done.'])
