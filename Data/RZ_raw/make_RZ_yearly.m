%load Ramey Zubairy news data and convert to yearly
% Note: The data used here are from Ramey and Zubairy (2018, JPE)
% "Government Spending Multipliers in Good Times and in Bad: Evidence from 
% US Historical Data". The file RZDAT.xlsx contains their raw data, and is
% available from their online appendix.


clear all
close all
clc


%load RZ quarterly news data
load RZ_data_raw
%news = nominal military news
%rgdp_pott6 = potential gdp measure
%pgdp = gdp deflator
%quarter = quarter

%nominal potential gdp
pngdp = rgdp_pott6.*pgdp;

%news as fraction of lagged nominal potential gdp
newsy = [NaN;news(2:end)./pngdp(1:end-1)];

%year
year = floor(quarter);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make yearly data

%year
year_y = unique(floor(quarter));
Ty = length(year_y);


news_y = zeros(Ty,1);
pngdp_y = zeros(Ty,1);
%loop over years
for ty = 1:Ty
    %sum of military announcements that year
    news_y(ty) = sum( news(year == year_y(ty)) );
    %sum of ngdp that year
    pngdp_y(ty) = sum( pngdp(year == year_y(ty)) );
end

%yearly news as fraction of ngdp
newsy_y = [NaN;news_y(2:end)./pngdp_y(1:end-1)];


plot(quarter,newsy)
hold on
plot(year_y,newsy_y)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export

%save, skipping first year since NaN
year_RZ = year_y(2:end);
newsy_RZ = newsy_y(2:end);

save RZ_yearly year_RZ newsy_RZ









