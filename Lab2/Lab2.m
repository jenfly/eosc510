% Lab 2 (14 Sep 2017): Excercise 2 w/ solutions 

% Units of variables given in the dataset:
%    elevation: meters above sea level
%    g: (mm water equivalent)/meter
%    continentality: degrees Celsius
%    cloud cover: percentage
%    temperature/ summer temperature: degrees Celsius
%    precipitation: mm/year
%    average winter precipitation: mm/month 

% NOTE: climate variables are given from climatological gridded data (CRU) provided for the whole globe on ca 50x50 km resolution
%
% Variables:
%                       g: [1x136 double]  glacier mass balance gradient (as an average over the observed period)
%        median_elevation: [1x136 double]  Median glacier elevation (derived from glacier hypsometry: glacier area vs elevation)
%      summer_temperature: [1x136 double]  Summer near-surface air temperature over a grid cell (ca 50x50 km) covering the glacier (averaged over the observed period)
%           precipitation: [1x136 double]  Total annaul precipitation over a grid cell (ca 50x50 km) covering the glacier (averaged over the observed period)
%                     lat: [1x136 double]  Geographical latitude of the glacier (in degrees)
%           max_elevation: [1x136 double]  Maximum glacier elevation
%                 WGMS ID: [1x136 int64]   Given glacier ID (from World Glacier Monitoring Service)
%    winter_precipitation: [1x136 double]  Same as for precipation but derived only over winter months 
%                     lon: [1x136 double]  Geographical longitude of the glacier (in degrees)
%             cloud_cover: [1x136 double]  Annual cloud cover over a grid cell (ca 50x50 km) covering the glacier (averaged over the observed period)
%          continentality: [1x136 double]  Maximum monthly temperature minus minimum monthly temperature over a year (over a grid cell covering the glacier and averaged over the observed period). The larger the continentality index the climate is more continental (rather than maritime). 

% LOAD THE DATA
load glaciers.mat;

% move the structure data into separate variables 
y=glaciers.lat;
x=glaciers.lon;
hmax=glaciers.max_elevation;
hmed=glaciers.median_elevation;
g=glaciers.g;
Pannual=glaciers.precipitation;
Pwinter=glaciers.winter_precipitation;
CI=glaciers.continentality;
cloud=glaciers.cloud_cover;
Tsummer=glaciers.summer_temperature;

% plot the world map with the locations of glaciers from the dataset
% if your Matlab version has the mapping toolbox use the following command:
% load coast 
% this loads lat and long of the continents coastlines globally
% otherwise use the file in your directory:
load coastline.mat

N=length(x);

figure;
plot(x,y,'ro'); 
hold on
plot(long,lat,'k');
title('glaciers with their g');

% getting rid of negative g (not theoretically possible)
[dummy index]=find(g > 0); % the indices of positive values are save in variable 'index' 

% place all the plausible predictor variables into a matrix 'xall' so that columns are variables and rows are glaciers
% transposing each variable b/c the original is given as a row but I want each as a column (vector) 
xall=[y(index)' x(index)' hmax(index)' hmed(index)' Pannual(index)' Pwinter(index)' CI(index)' cloud(index)' Tsummer(index)'];

% place mass balance gradianet variable into a vector 'Yg'
Yg=g(index)';


