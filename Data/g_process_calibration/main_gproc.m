%Script to calibrate the g process using an AR(1) regression on linearly
%detrended log real government spending.

clear all; 
clc; 
close all;

table=table2array(readtable('FREDRealGovernmentConsumptionExpenditures.csv'));
g      =table(:,2);

years = 1971:2013;

%detrend with linear trend
log_g_stat = detrend(log(g));
%AR(1) regression
log_g_Y = log_g_stat(2:end);
log_g_X = log_g_stat(1:end-1);
[rho_g, ~, res_g] = regress(log_g_Y, log_g_X);

%parameters used in calibration
rho_g
sigma_inn_g = std(res_g)






