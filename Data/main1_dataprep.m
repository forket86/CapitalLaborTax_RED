%Script to load and analyze annual tax, spending, gdp, and debt data
%This script loads the raw data which has been manually collected and saved
%into data.mat. It processes the data and produces the following outputs:
% 1. produces data_out.csv for use in local projection estimation file
% main2_estimateLP.do
% 2. produces Figure D1 and associated numbers
% 3. produces data for model estimation and comparison used in final
% columns of Tables 2 and 5


clear all
close all
clc


%set = 1 to save figures in /Figures folder, or = 0 not to save
save_figs = 1
figsSavePath = '../Figures/';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data and extract variables

%load data (raw NIPA tables, row # corresponds to row number on NIPA table and website)
load data.mat

% number of time periods in NIPA dataset
T = length(date);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract variables
% NIPA: millions of dollars
% Debt: billions of dollars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLES 1.1.4 - 1.1.5

% Output price deflator
DEF = table1p1p4(1,:)';

% Nominal GDP
NGDP = table1p1p5(1,:)';

% Nominal consumption
PCE = table1p1p5(2,:)';

% Nominal government spending and investment
GSP = table1p1p5(22,:)';

% Nominal government defense spending and investment
DFN = table1p1p5(24,:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLE 1.12

% Nominal compensation of employees
CEM = table1p12(2,:)';

% Nominal wages and salary accruals
WSA = table1p12(3,:)';

% Nominal proprietors income
PRI = table1p12(9,:)';

% Nominal rental income
RI = table1p12(12,:)';

% Nominal corporate profits
CP = table1p12(13,:)';

% Nominal interest income
NI = table1p12(18,:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TABLES 3.1-3.3

% Nominal taxes on production and imports
TPI = table3p1(4,:)';

% Nominal taxes on corporate income
CT = table3p1(5,:)';

% Nominal contributions to social security
CSI = table3p1(7,:)';

% Nominal personal income tax (federal, state and local)
PIT = table3p2(3,:)' + table3p3(4,:)';

% Nominal state and local property taxes
PRT = table3p3(9,:)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ramey Zubairy news shock data

load RZ_yearly.mat

%merge data and add NaN for missing ranges
newsy = NaN(T,1);
for t = 1:T
    today_ = date(t);
    t1d = find(year_RZ==today_);
    if ~isempty(t1d)
        newsy(t) = newsy_RZ(t1d);
    end
end
clear newsy_RZ year_RZ


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debt data

%convert debt to millions of dollars
debt = debt*1000;
debt_public = debt_public*1000;
debt_public_cbo = debt_public_cbo*1000;

%merge debt data and add NaN for missing range (FRED DEBT DATA)
debt_ = NaN(T,1);
debt_public_ = NaN(T,1);
for t = 1:T
    today_ = date(t);
    t1d = find(date_debt==today_);
    if ~isempty(t1d)
        debt_(t) = debt(t1d);
        debt_public_(t) = debt_public(t1d);
    end
end
debt = debt_;
debt_public = debt_public_;
clear debt_ debt_public_ date_debt


%merge debt data and add NaN for missing range (CBO DEBT DATA)
debt_public_cbo_ = NaN(T,1);
for t = 1:T
    today_ = date(t);
    t1d = find(date_cbo==today_);
    if ~isempty(t1d)
        debt_public_cbo_(t) = debt_public_cbo(t1d);
    end
end
debt_public_cbo = debt_public_cbo_;
clear debt_public_cbo_ date_cbo


%splice CBO debt data onto beginning of FRED debt data series, scaling so
%that values in first year of FRED data are the same
td1 = find(~isnan(debt_public),1,'first');
debt_public(1:td1) = debt_public_cbo(1:td1)*debt_public(td1)/debt_public_cbo(td1);




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD OUR VARIABLES


% Consumption tax rate
tauC = (TPI-PRT)./(PCE-(TPI-PRT));

% Labour and capital personal income
LI = WSA + PRI/2;
CI = PRI/2 + RI + CP + NI;

% Personal income tax rate
tauP = PIT./(LI + CI);

% Labour tax rate
tauL = (tauP.*LI + CSI)./(CEM + PRI/2);

% Capital tax rate
tauK = (tauP.*CI + CT + PRT)./(CI + PRT);

% Real GDP
y = NGDP./DEF;

% Real government spending
g = GSP./DEF;
g_dfn = DFN./DEF;

% real gov debt
rdebt = debt./DEF;
rdebt_pub = debt_public./DEF;
%rdebt_pub_cbo = debt_public_cbo./DEF;

% debt to GDP
debty = debt./NGDP;
debty_pub = debt_public./NGDP;
%debty_pub_cbo = debt_public_cbo./NGDP;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BIG PLOT OF ALL VARIABLES


%time vector for plots, with yearly data put in midpoint (June) of the year
dates = datetime(date,6,30);


%optional: big plot to look at all variables
make_bigPlot = 0;
if make_bigPlot

    figsize = [100,100,700,500]
    figure('position',figsize)


    MM = 3; NN = 3; ip = 1;
    xmin = min(date);
    xmax = max(date);
    xmin = datetime(xmin,6,30);
    xmax = datetime(xmax,6,30);

    subplot(MM,NN,ip)
    plot(dates,tauK,'linewidth',2)
    title('Capital Tax \tau_K')
    grid on
    ylims = ylim;
    recessionplot
    ylim(ylims)
    xlim([xmin,xmax])
    ip = ip + 1;

    subplot(MM,NN,ip)
    plot(dates,tauL,'linewidth',2)
    title('Labor Tax \tau_L')
    grid on
    ylims = ylim;
    recessionplot
    ylim(ylims)
    xlim([xmin,xmax])
    ip = ip + 1;

    subplot(MM,NN,ip)
    plot(dates,tauP,'linewidth',2)
    title('Private Income Tax \tau_P')
    grid on
    ylims = ylim;
    recessionplot
    ylim(ylims)
    xlim([xmin,xmax])
    ip = ip + 1;

    subplot(MM,NN,ip)
    plot(dates,tauC,'linewidth',2)
    title('Consumption Tax \tau_C')
    grid on
    ylims = ylim;
    recessionplot
    ylim(ylims)
    xlim([xmin,xmax])
    ip = ip + 1;

    subplot(MM,NN,ip)
    plot(dates,g,'linewidth',2)
    title('Gov Spending g')
    grid on
    ylims = ylim;
    recessionplot
    ylim(ylims)
    xlim([xmin,xmax])
    ip = ip + 1;

    subplot(MM,NN,ip)
    plot(dates,debty,'linewidth',2)
    hold on
    plot(dates,debty_pub,'linewidth',2)
    title('Debt/GDP')
    grid on
    ylims = ylim;
    xlim([xmin,xmax])
    recessionplot
    ylim(ylims)
    legend('Debt','Debt held by public')
    xlim([xmin,xmax])
    ip = ip + 1;

    subplot(MM,NN,ip)
    plot(dates,g./y,'linewidth',2)
    title('Gov Spending/GDP')
    grid on
    ylims = ylim;
    recessionplot
    ylim(ylims)
    xlim([xmin,xmax])
    ip = ip + 1;

    subplot(MM,NN,ip)
    plot(dates,y,'linewidth',2)
    title('Real GDP y')
    grid on
    ylims = ylim;
    recessionplot
    ylim(ylims)
    xlim([xmin,xmax])
    ip = ip + 1;

    % subplot(MM,NN,ip)
    % plot(dates,g_dfn,'linewidth',2)
    % title('Real Gov Def Spending')
    % grid on
    % ylims = ylim;
    % recessionplot
    % ylim(ylims)
    % xlim([xmin,xmax])
    % ip = ip + 1;

    subplot(MM,NN,ip)
    plot(dates,newsy,'linewidth',2)
    title('RZ defense news shock')
    grid on
    ylims = ylim;
    recessionplot
    ylim(ylims)
    xlim([xmin,xmax])
    ip = ip + 1;

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE D1 AND ASSOCIATED NUMBERS (COMPOSITION OF TAUL AND TAUK)


% Personal income tax rate
mtauP = mean(tauP);

% Labour tax rate
tauL_cf = (mtauP.*LI + CSI)./(CEM + PRI/2);

% Capital tax rate
tauK_cf = (mtauP.*CI + CT + PRT)./(CI + PRT);


%time vector for plots, with yearly data put in midpoint (June) of the year
dates = datetime(date,6,30);


figsize = [100,100,600,200]
figure('position',figsize)


MM = 1; NN = 2; ip = 1;
xmin = min(date);
xmax = max(date);
xmin = datetime(xmin,6,30);
xmax = datetime(xmax,6,30);

subplot(MM,NN,ip)
plot(dates,tauK,'linewidth',2)
hold on
plot(dates,tauK_cf,'linewidth',2)
title('Capital Tax \tau_K')
grid on
ylims = ylim;
recessionplot
ylim(ylims)
xlim([xmin,xmax])
legend('true','fix \tau_p','location','southeast')
ip = ip + 1;

subplot(MM,NN,ip)
plot(dates,tauL,'linewidth',2)
hold on
plot(dates,tauL_cf,'linewidth',2)
title('Labor Tax \tau_L')
grid on
ylims = ylim;
recessionplot
ylim(ylims)
xlim([xmin,xmax])
legend('true','fix \tau_p','location','southeast')
ip = ip + 1;


if save_figs
    savefig('../Figures/FigureD1.fig')
end




% Average contributions

% Labour tax rate
av_tPLI = mean(tauP.*LI);
av_CSI = mean(CSI);
tauLfrac_tPLI = av_tPLI/(av_tPLI+av_CSI)
tauLfrac_CSI = av_CSI/(av_tPLI+av_CSI)

% Capital tax rate
av_tPLC = mean(tauP.*CI);
av_CT = mean(CT);
av_PRT = mean(PRT);
tauKfrac_tPLC = av_tPLC/(av_tPLC+av_CT+av_PRT)
tauKfrac_CT = av_CT/(av_tPLC+av_CT+av_PRT)
tauKfrac_PRT = av_PRT/(av_tPLC+av_CT+av_PRT)





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA EXPORT FOR LOCAL PROJECTIONS IN STATA


%create empty table
Tab = table;

%load variables into table
varlist = {'date','y','g','tauK','tauL','tauP','tauC','rdebt','rdebt_pub','newsy'};
for i = 1:length(varlist)
    eval(['Tab.',varlist{i},' = ',varlist{i}])
end

%delete 2016 onwards since no news shock data
Tab(Tab.date>=2016,:) = [];

%export to csv
writetable(Tab,'data_out.csv')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BUILD DATA FOR CALIBRATION TABLES: TABLE 2 AND 5

%data range
mindate1 = 1971;
maxdate1 = 2013;
ind = (date >= mindate1) & (date <= maxdate1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AVERAGES

%g/y
gy = g./y;

%loop over time ranges
E_gy = mean(gy(ind));
E_tauK = mean(tauK(ind));
E_tauL = mean(tauL(ind));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STDS AND AUTOCORRS

% Detrending method: do not detrend since tauK and tauL are rates which
% have no trend in the model
dtm = 'off';

% log(1+tauL)
ltauL = detrender(log(1+tauL(ind)),dtm);

% log(1+tauK)
ltauK = detrender(log(1+tauK(ind)),dtm);

% tauK
data = ltauK;
std_tauK = std(data);
rho_tauK = corr(data(2:end),data(1:end-1));

% tauL
data = ltauL;
std_tauL = std(data);
rho_tauL = corr(data(2:end),data(1:end-1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE FINAL COLUMN OF TABLES 2 AND 5


FinalCol_Tab2and5 = [E_tauK;E_tauL;std_tauK;std_tauL;rho_tauK;rho_tauL]














