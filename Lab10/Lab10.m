%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 10 (9 Nov 2016)
% MLP NN modeling on a synthetic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
load coast % this loads lat and long of the continents coastlines 

figure;
plot(x,y,'ro'); 
hold on
plot(long,lat,'k');
title('glaciers with their g');

% getting rid of negative g (not theoretically possible)
[dummy index]=find(g > 0); % the indices of positive values are save in variable 'index' 

% place all the variables into a matrix 'xall' so that columns are varibales and rows are glaciers
% also: take absolute value of y (latitudes) since we are interested only in how far north or south is the glacier

Yg=g(index)';
xall=[abs(y(index))' x(index)' hmax(index)' hmed(index)' Pannual(index)' Pwinter(index)' CI(index)' cloud(index)' Tsummer(index)'];

% apply PCA on the data with 9 variables and n glaciers

% standardizing the variables (x1-mean(x1))/std(x1), etc
mean_xall=mean(xall);
std_xall=std(xall);

n=length(index);
xall_mean=repmat(mean_xall,n,1);
xall_std=repmat(std_xall,n,1);

xalls=(xall-xall_mean)./xall_std;

[eigenvectors,PCs,eigval]=princomp(xalls);

figure;
subplot(1,2,1)
plot(xalls(1,:));

subplot(1,2,2)
plot(xalls(2,:));

% contribution of each mode to total variance:
variance=eigval./sum(eigval);

figure; 
plot([1:length(variance)],variance,'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');
% from this figure I decided to keep 3 first modes

% plot the first 3 modes (eigenvectors and PCs)
figure; 
for i=1:3
subplot(3,2,2*i-1)
plot(eigenvectors(:,i));
title(['e_',num2str(i)]);
xlabel('variables');

subplot(3,2,2*i)
plot(PCs(:,i));
xlabel('index of glaciers');
title(['PC',num2str(i)]);
end





