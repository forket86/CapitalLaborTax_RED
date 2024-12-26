%Script to plot local projection results from Stata.
%This script loads the output files with name structure 
%lp_results_sample*_lags*.csv which are produced by main2_estimateLP.do in
%Stata. This script produces:
% 1. Figure 6: main LP reults
% 2. Figure D2: LP output result


clear all
close all
clc


%set = 1 to save figures in /Figures folder, or = 0 not to save
save_figs = 1
figsSavePath = '../Figures/';


%sample choice. 1 = 1929 to 2015. 3 = 1947 to 1995
%Baseline used in paper is sample = 1, lags = 2
sample = 1;
%lag choice. 1 or 2
lags = 2;
%scale size of shock so that peak gov spending response is around 50%
scale = .221;

%load the data for the chosen sample and lags
data = readtable(['lp_results_sample',num2str(sample),'_lags',num2str(lags),'.csv']);

%extract table into vectors
for iv = 1:length(data.Properties.VariableNames)
    %extract 50-65 data, story as X_long
    eval([data.Properties.VariableNames{iv},' = data{:,iv};'])
end
clear data data_ws iv





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 6

sig = 1.96;
lcol = [0 0.4470 0.7410]; %main line colour
fcol = 0.9*[1,1,1]; %std range fill colour
falpha = 0.5; %std range transparency
flcol = 0.7*[1,1,1]; %std range edge colour

xlims = [0,4];
xtks = (min(xlims):max(xlims));

figsize = [100,100,700,200];

MM = 2; NN = 2; ip = 1;
%figure('position',figsize)
figure

subplot(MM,NN,ip)
%%% Set variable here
Y = scale*100*b_lg;
Ystd = scale*100*se_lg;
ttle = 'Real Gov Spending';
%%% Don't edit here
Yup = Y + sig*Ystd;
Ylow = Y - sig*Ystd;
fill([h; flipud(h)],[Yup; flipud(Ylow)],fcol,'FaceAlpha',falpha,'EdgeColor', flcol)
hold on
plot(h,Y,'linewidth',2,'color',lcol)
grid on
box off
%title(ttle)
xlabel('years')
xlim(xlims)
xticks(xtks)
ylabel(['\Delta ',ttle,', %'])
%%%
ip = ip + 1;


subplot(MM,NN,ip)
%%% Set variable here
Y = scale*100*b_tauk;
Ystd = scale*100*se_tauk;
ttle = 'Capital tax';
%%% Don't edit here
Yup = Y + sig*Ystd;
Ylow = Y - sig*Ystd;
fill([h; flipud(h)],[Yup; flipud(Ylow)],fcol,'FaceAlpha',falpha,'EdgeColor', flcol)
hold on
plot(h,Y,'linewidth',2,'color',lcol)
grid on
box off
%title(ttle)
xlabel('years')
xlim(xlims)
xticks(xtks)
ylabel(['\Delta ',ttle,', p.p.'])
%%%
ip = ip + 1;


subplot(MM,NN,ip)
%%% Set variable here
Y = scale*100*b_taul;
Ystd = scale*100*se_taul;
ttle = 'Labor tax';
%%% Don't edit here
Yup = Y + sig*Ystd;
Ylow = Y - sig*Ystd;
fill([h; flipud(h)],[Yup; flipud(Ylow)],fcol,'FaceAlpha',falpha,'EdgeColor', flcol)
hold on
plot(h,Y,'linewidth',2,'color',lcol)
grid on
box off
%title(ttle)
xlabel('years')
xlim(xlims)
xticks(xtks)
ylabel(['\Delta ',ttle,', p.p.'])
%%%
ip = ip + 1;


subplot(MM,NN,ip)
%%% Set variable here
Y = scale*100*b_lb;
Ystd = scale*100*se_lb;
ttle = 'Real Debt';
%%% Don't edit here
Yup = Y + sig*Ystd;
Ylow = Y - sig*Ystd;
fill([h; flipud(h)],[Yup; flipud(Ylow)],fcol,'FaceAlpha',falpha,'EdgeColor', flcol)
hold on
plot(h,Y,'linewidth',2,'color',lcol)
grid on
box off
%title(ttle)
xlabel('years')
xlim(xlims)
xticks(xtks)
ylabel('%')
ylabel(['\Delta ',ttle,', %'])
%%%
ip = ip + 1;

if save_figs
    savefig('../Figures/Figure6.fig')
end





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE D2



figsize = [100,100,200,150];
figure('position',figsize)
%%% Set variable here
Y = scale*100*b_ly;
Ystd = scale*100*se_ly;
ttle = 'Real GDP';
%%% Don't edit here
Yup = Y + sig*Ystd;
Ylow = Y - sig*Ystd;
fill([h; flipud(h)],[Yup; flipud(Ylow)],fcol,'FaceAlpha',falpha,'EdgeColor', flcol)
hold on
plot(h,Y,'linewidth',2,'color',lcol)
grid on
box off
%title(ttle)
xlabel('years')
xlim(xlims)
xticks(xtks)
ylabel(['\Delta ',ttle,', %'])
%%%
ip = ip + 1;



if save_figs
    savefig('../Figures/FigureD2.fig')
end




