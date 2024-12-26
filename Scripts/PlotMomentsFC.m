
%% Fig 3 e 4 comparison with Fahri
clear all;
clc;
close all;

rng(1)

%% SET THIS TO 1 TO RECALCULATE MOMENTS
RECALCULATE_SIMULATION=0;
infinityHorizon=5000;

if RECALCULATE_SIMULATION
    
    addpath('../Utils/');
    addpath('./Measurability Neural Integrand/');
    
    cutInit=50;
    costgrid=[0.1:0.1:0.9 5 20 30 39 50];
    for costIndex=1:length(costgrid)
        
        if costgrid(costIndex)<1
            load(['./Measurability Neural Integrand/Good/Curve/Both/Curve0dot' num2str(costIndex)])
        else
            load(['./Measurability Neural Integrand/Good/Curve/Both/Curve' num2str(costgrid(costIndex))])
        end
        
        T=infinityHorizon;
        rng(1)
        %storetauK=tauK;
        Simulate
        %plot([storetauK(1:t)' tauK'])
        
        ETauL(costIndex)=mean(tauL(cutInit:T));
        ETauK(costIndex)=mean(tauK(cutInit:T));
        StdTauL(costIndex)=std(log(tauL(cutInit:T)+1));
        StdTauK(costIndex)=std(log(tauK(cutInit:T)+1));
        
        if costgrid(costIndex)<1
            load(['./Measurability Neural Integrand/Good/Curve/Lowchi/CurveLOWCHI0dot' num2str(costIndex)])
        else
            load(['./Measurability Neural Integrand/Good/Curve/Lowchi/CurveLOWCHI' num2str(costgrid(costIndex))])
        end
        
        T=infinityHorizon;
        rng(1)
        %storetauK=tauK;
        Simulate
        %plot([storetauK(1:t)' tauK'])
        
        ETauL2(costIndex)=mean(tauL(cutInit:T));
        ETauK2(costIndex)=mean(tauK(cutInit:T));
        StdTauL2(costIndex)=std(log(tauL(cutInit:T)+1));
        StdTauK2(costIndex)=std(log(tauK(cutInit:T)+1));
    end
else
    load('Simulations')
end

%% figure
close all
hline_base = 3;
hline_data = 1;
hline_comp = 2;
fsize = 14;
if 0
    figure
    subplot(2,2,1)
    hold on
    plot([0.01 costgrid],[mean(tauK_FC(cutInit:end)) ETauK], 'linewidth',hline_base)
    plot([0.01 costgrid],[mean(tauK_FC(cutInit:end)) ETauK2],'--','linewidth',hline_base)
    plot([0.01 costgrid],0.3665*ones(1,length(costgrid)+1),'k-.','linewidth',hline_data)
    title('$E(\tau^k)$', 'fontsize', fsize,'interpreter','latex')
    xlabel('$\gamma$', 'fontsize', fsize,'interpreter','latex')
    grid on
    xlim([0,50])
    ylim([-2*10^-3,0.4])
    
    subplot(2,2,2)
    
    hold on
    plot([0.01 costgrid],[mean(tauL_FC(cutInit:end)) ETauL],'linewidth',hline_base)
    plot([0.01 costgrid],[mean(tauL_FC(cutInit:end)) ETauL2],'--','linewidth',hline_base)
    plot([0.01 costgrid],0.2250*ones(1,length(costgrid)+1),'k-.','linewidth',hline_data)
    title('$E(\tau^l)$', 'fontsize', fsize,'interpreter','latex')
    xlabel('$\gamma$', 'fontsize', fsize,'interpreter','latex')
    grid on
    xlim([0,50])
    ylim([0,0.4])
end
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
ylim([0,0.07])
