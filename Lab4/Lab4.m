% PCA on glacier data (28 Sep 2017)
% Units for variables:
%    elevation: meters above sea level
%    g: (mm water equivalent)/meter
%    continentality: degrees Celsius
%    cloud cover: percentage
%    temperature/ summer temperature: degrees Celsius
%    precipitation: mm/year
%    average winter precipitation: mm/month 
%
% Variables:
%                       g: [1x136 double]  glacier mass balance gradient (as an averege over observed period)
%        median_elevation: [1x136 double]  Median glacier elevation
%      summer_temperature: [1x136 double]  Summer near-surface air temperature over a grid cell (ca 50x50 km) covering the glacier (averaged over the observed period)
%           precipitation: [1x136 double]  Total annaul precipitation over a grid cell (ca 50x50 km) covering the glacier (averaged over the observed period)
%                     lat: [1x136 double]  Geographical latitude of the glacier (in degrees)
%           max_elevation: [1x136 double]  Maximum glacier elevation
%                 WGMS ID: [1x136 int64]   Given glacier ID (from World Glacier Monitoring Service)
%    winter_precipitation: [1x136 double]  Same as for precipation but derived only over winter months 
%                     lon: [1x136 double]  Geographical longitude of the glacier (in degrees)
%             cloud_cover: [1x136 double]  Annual cloud cover over a grid cell (ca 50x50 km) covering the glacier (averaged over the observed period)
%          continentality: [1x136 double]  Maximum monthly temperature minus minimum monthly temperature over a year over a grid cell covering the glacier (averaged over the observed period) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES:
% 1) Use all the variable above except WGMS ID -> done below in the code 
% 2) Remove the glaciers with negative g (not theoretically possible) -> done below in the code
% 3) The variables have different scales so some pre-processing (hint: standardization) is recommended. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% load coast 
% this loads lat and long of the continents coastlines globally
% otherwise use the file in your directory:
load coastline.mat

N=length(x);
figure;
plot(x,y,'ro'); 
hold on
plot(long,lat,'k');
title('glaciers in the sample');

% getting rid of negative g (not theoretically possible)
[dummy index]=find(g > 0); % the indices of positive values are save in variable 'index' 

% place all the variables into a matrix 'xall' so that columns are variables and rows are glaciers
% here we removed the variable WGMI ID since this one is not needed in the analysis
xall=[g(index)' y(index)' x(index)' hmax(index)' hmed(index)' Pannual(index)' Pwinter(index)' CI(index)' cloud(index)' Tsummer(index)'];


