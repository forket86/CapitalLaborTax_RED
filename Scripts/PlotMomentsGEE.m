
%% Fig 3 e 4 comparison with Fahri
clear all;
clc;
close all;

rng(1)

%% SET THIS TO 1 TO RECALCULATE MOMENTS
RECALCULATE_SIMULATION=1;
infinityHorizon=1000;
figure 
hold on
if RECALCULATE_SIMULATION
    
    addpath('../Utils/');
    addpath('./Measurability Neural Integrand/');
    
    cutInit=50;
    costgrid=[2 6 39];
    for costIndex=1:length(costgrid)
        
        if costgrid(costIndex)<1
            load(['./Measurability Neural Integrand/Good/CurveGEE/Both/Curve0dot' num2str(10*costgrid(costIndex))])
        else
            load(['./Measurability Neural Integrand/Good/CurveGEE/Both/Curve' num2str(costgrid(costIndex))])
        end
        
        T=infinityHorizon;
        rng(1)
        %storetauK=tauK;
        %Simulate
        %plot([storetauK(1:t)' tauK'])
        subplot(2,1,1)
        plot(g_shocks(cutInit:T))
        subplot(2,1,2)
        hold on
        plot(tauK(cutInit:T))
        ETauL(costIndex)=mean(tauL(cutInit:T));
        ETauK(costIndex)=mean(tauK(cutInit:T));
        StdTauL(costIndex)=std(log(tauL(cutInit:T)+1));
        StdTauK(costIndex)=std(log(tauK(cutInit:T)+1));
        
        Y=Z*k(1:T).^alpha.*labor(1:T).^(1-alpha);
        wage=(1-alpha)*Y(1:T);
        r=alpha*Y(1:T)-delta;
        StdTauLIncome(costIndex)=std(log(labor(cutInit:T).*wage(cutInit:T).*tauL(cutInit:T)));
        StdTauKIncome(costIndex)=std(log(k(cutInit:T).*r(cutInit:T).*tauK(cutInit:T)));
        
        Perror(1,costIndex)=rsquare(Output1(1,start:T)',Poutput(1,:)');
        Perror(2,costIndex)=rsquare(Output2(1,start:T)',Poutput(2,:)');
        Perror(3,costIndex)=rsquare(Output3(1,start:T)',Poutput(3,:)');
        Perror(4,costIndex)=rsquare(Output4(1,start:T)',Poutput(4,:)');
        Perror(5,costIndex)=rsquare(Output5(1,start:T)',Poutput(5,:)');
        
        
        
%         if costgrid(costIndex)<1
%             load(['./Measurability Neural Integrand/Good/Curve/Lowchi/CurveLOWCHI0dot' num2str(costIndex)])
%         else
%             load(['./Measurability Neural Integrand/Good/Curve/Lowchi/CurveLOWCHI' num2str(costgrid(costIndex))])
%         end
%         
%         T=infinityHorizon;
%         rng(1)
%         %storetauK=tauK;
%         Simulate
%         %plot([storetauK(1:t)' tauK'])
%         
%         ETauL2(costIndex)=mean(tauL(cutInit:T));
%         ETauK2(costIndex)=mean(tauK(cutInit:T));
%         StdTauL2(costIndex)=std(log(tauL(cutInit:T)+1));
%         StdTauK2(costIndex)=std(log(tauK(cutInit:T)+1));
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

figure
subplot(2,2,1)
hold on
plot(costgrid,ETauK, 'linewidth',hline_base)
%plot(costgrid,ETauK2,'--','linewidth',hline_base)
plot(costgrid,0.3665*ones(1,length(costgrid)),'k-.','linewidth',hline_data)
title('$E(\tau^k)$', 'fontsize', fsize,'interpreter','latex')
xlabel('$\gamma$', 'fontsize', fsize,'interpreter','latex')
grid on
xlim([min(costgrid),max(costgrid)])
ylim([min([ETauK ETauL]),max([ETauK ETauL 0.3665])])

subplot(2,2,2)

hold on
plot(costgrid,ETauL,'linewidth',hline_base)
%plot(costgrid,ETauL2,'--','linewidth',hline_base)
plot(costgrid,0.2250*ones(1,length(costgrid)),'k-.','linewidth',hline_data)
title('$E(\tau^l)$', 'fontsize', fsize,'interpreter','latex')
xlabel('$\gamma$', 'fontsize', fsize,'interpreter','latex')
grid on
xlim([min(costgrid),max(costgrid)])
ylim([min([ETauK ETauL]),max([ETauK ETauL])])

subplot(2,2,3)

hold on
plot(costgrid,StdTauK,'linewidth',hline_base)
%plot(costgrid,StdTauK2,'--','linewidth',hline_base)
plot(costgrid,0.0212*ones(1,length(costgrid)),'k-.','linewidth',hline_data)
title('$\sigma(\tau^k)$', 'fontsize', fsize,'interpreter','latex')
xlabel('$\gamma$', 'fontsize', fsize,'interpreter','latex')
grid on
xlim([min(costgrid),max(costgrid)])
ylim([min([StdTauK StdTauL]),max([StdTauK StdTauL])])

subplot(2,2,4)

hold on
plot(costgrid,StdTauL,'linewidth',hline_base)
%plot(costgrid,StdTauL2,'--','linewidth',hline_base)
plot(costgrid,0.0159*ones(1,length(costgrid)),'k-.','linewidth',hline_data)
title('$\sigma(\tau^l)$', 'fontsize', fsize,'interpreter','latex')
xlabel('$\gamma$', 'fontsize', fsize,'interpreter','latex')
grid on
xlim([min(costgrid),max(costgrid)])
ylim([min([StdTauK StdTauL]),max([StdTauK StdTauL])])

figure
subplot(5,1,1)
plot(costgrid,Perror(1,:))
grid on
xlabel('$\gamma$', 'fontsize', fsize-2,'interpreter','latex')
%title('$R^2$ prediction 1', 'fontsize', fsize,'interpreter','latex')
subplot(5,1,2)
plot(costgrid,Perror(2,:))
xlabel('$\gamma$', 'fontsize', fsize-2,'interpreter','latex')
%title('$R^2$ prediction 2', 'fontsize', fsize,'interpreter','latex')
grid on
subplot(5,1,3)
plot(costgrid,Perror(3,:))
xlabel('$\gamma$', 'fontsize', fsize-2,'interpreter','latex')
%title('$R^2$ prediction 3', 'fontsize', fsize,'interpreter','latex')
grid on
subplot(5,1,4)
plot(costgrid,Perror(4,:))
xlabel('$\gamma$', 'fontsize', fsize-2,'interpreter','latex')
%title('$R^2$ prediction 4', 'fontsize', fsize,'interpreter','latex')
grid on
subplot(5,1,5)
plot(costgrid,Perror(5,:))
grid on
xlabel('$\gamma$', 'fontsize', fsize-2,'interpreter','latex')
%title('$R^2$ prediction 5', 'fontsize', fsize,'interpreter','latex')