% PCA on glacier data (28 Sep 2017) with solutions 

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
load coast % this loads lat and long of the continents coastlines 

N=length(x);
figure;
plot(x,y,'ro'); 
hold on
plot(long,lat,'k');
title('glaciers with their g');

% getting rid of negative g (not theoretically possible)
[dummy index]=find(g > 0); % the indices of positive values are save in variable 'index' 

% place all the variables into a matrix 'xall' so that columns are varibales and rows are glaciers
xall=[g(index)' y(index)' x(index)' hmax(index)' hmed(index)' Pannual(index)' Pwinter(index)' CI(index)' cloud(index)' Tsummer(index)'];

% 1) What are the most characteristic features in the space of these variables across the glaciers?
% to answer this we apply PCA on the data with 10 variables and 136 glaciers
% the eignevectors (which carry most of the variance) will represent these characteristic features 

% standardizing the variables (x1-mean(x1))/std(x1), etc
mean_xall=mean(xall);
std_xall=std(xall);

n=length(index);
xall_mean=repmat(mean_xall,n,1);
xall_std=repmat(std_xall,n,1);

xalls=(xall-xall_mean)./xall_std;

[eigenvectors,PCs,eigval]=princomp(xalls);

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
subplot(2,3,i)
plot(eigenvectors(:,i));
title(['e_',num2str(i)]);
xlabel('variables');

subplot(2,3,i+3)
plot(PCs(:,i));
xlabel('index of glaciers');
title(['PC',num2str(i)]);
end

% 2) Can the glaciers be characterized with fewer variables (modes)?
% Yes, we can reconstruct the data using, for example, the first 2 to 4 modes 
% Below I choose the first 3 modes:
% reconstructed data from the first 3 modes
invV=inv(eigenvectors');
x_rec=(invV(:,1:3)*PCs(:,1:3)')';

% 3) If yes, how do the most 'generic' glaciers look like according to these fewer modes?
% finding the 'generic' glacier, i.e. the one with minimum RMSE between reconstructed and original data
err=xalls-x_rec;
RMSE=(sum((err.^2)')./100).^0.5;
[minRMSE indmin]=sort(RMSE,'ascend');

% plot the most generic glacier characteristics
figure;
plot(xalls(indmin(1),:));
title('generic glacier characteristics');
xlim([0 11])
set(gca,'XTick',[1:1:10],'XTicklabel',{'g','lat','lon','hmax','hmed','P_{ann}','P{winter}','CI','cloud','T_{summer}'})

% highlight five most generic glaciers on the map
figure;
plot(x,y,'bo',xall(indmin(1:5),3),xall(indmin(1:5),2),'r*');
legend('all glaciers','generic glaciers')
hold on
plot(long,lat,'k');

% 4) What happens when you transpose the input matrix to PCA? Do the original eignevectors become
% PCs, and the original PCs become eigenvectors?

% doing PCA on trensposed inout data
[eigenvectors,PCs,eigval]=princomp(xalls');

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
subplot(2,3,i)
plot(eigenvectors(:,i));
title(['e_',num2str(i)]);
xlabel('variables');

subplot(2,3,i+3)
plot(PCs(:,i));
xlabel('index of glaciers');
title(['PC',num2str(i)]);
end

% comparing this figure with the one derived from the original data we can see that,
% apart from the mode 1, the second two PCs resemble the eigenvectors from the original data


