% Lab 2 (14 Sep 2017): Excercise with solutions 

% Units of variables given in the dataset:
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

% place all the predictor variables into a matrix 'xall' so that columns are varaibles and rows are glaciers
xall=[y(index)' x(index)' hmax(index)' hmed(index)' Pannual(index)' Pwinter(index)' CI(index)' cloud(index)' Tsummer(index)'];

% place g variable into a vector 'Yg'
Yg=g(index)';

% SOLUTIONS:

% 1) apply MLR and stepwise regression on standardized input
mean_xall=mean(xall);
std_xall=std(xall);

% standardizing the predictors (x1-mean(x1))/std(x1), etc
n=size(Yg,1);
xall_mean=repmat(mean_xall,n,1);
xall_std=repmat(std_xall,n,1);

xall_standard=(xall-xall_mean)./xall_std;

%regression with standardized variables
X=[ones(n,1) xall_standard];   
[b2,BINT2,R2,RINT2,STATS2]=regress(Yg,X);
y_regrMLR=X*b2;

% stepwise regression with standardized values
[a2 SE PVAL INMODEL STATS NEXTSTEP HISTORY]=stepwisefit(xall_standard,Yg,'penter',0.05);

% regression coefficients
a2=a2

% constant coefficient 
a20=STATS.intercept

% take out only the relevant predictors (from stepwisefit)
[ind1 ind2]=find(INMODEL == 1);
b=regress(Yg,X(:,[1 ind2+1]))
y_regrSTEP=X(:,[1 ind2+1])*b;

[r p]=corrcoef([Yg y_regrMLR y_regrSTEP]);

xline=[min([Yg;y_regrMLR]):max([Yg;y_regrMLR])]; 
figure;
subplot(2,1,1)
plot(Yg,y_regrMLR,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Yg;y_regrMLR]) max([Yg;y_regrMLR])]);
ylim([min([Yg;y_regrMLR]) max([Yg;y_regrMLR])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['MLR r= ' num2str(r(1,2),'%2.2f')]);

subplot(2,1,2)
plot(Yg,y_regrSTEP,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Yg;y_regrMLR]) max([Yg;y_regrMLR])]);
ylim([min([Yg;y_regrMLR]) max([Yg;y_regrMLR])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['Stepwise regression r= ' num2str(r(1,3),'%2.2f')]);

% Option 2 (calibration and validation)
% Calibrate and validate
% set an arrey from 1 to the total number of predictors (here = 4)
v0 = [1:size(xall,2)];

for kk=1:size(xall,2)
Ch = nchoosek(v0,kk);

for j=1:length(Ch(:,1))
% calibrate:
B = regress(Yg(1:2:end,:),X(1:2:end,[1 Ch(j,:)+1]));
% validate:
Ytest=X(2:2:end,[1 Ch(j,:)+1])*B;
[rtest0 ptest0]=corrcoef([Yg(2:2:end,:) Ytest]);
rtest(j)=rtest0(1,2);
ptest(j)=ptest0(1,2);
end
[pbest(kk) index(kk)]=min(ptest);
display('combination with columns:')
Ch(index(kk),:)
C(kk).predictors=Ch(index(kk),:);
display('r=') 
Rbest(kk)=rtest(index(kk))
clear ptest
clear rtest
end

% find the best combination and plot the final model
[pbestfinal indexfinal]=min(pbest);  % the best model is the one with the smallest p value
% recalculate the coefficents for this model:
Bfinal = regress(Yg(1:2:end,:),X(1:2:end,[1 C(indexfinal).predictors+1]));
% derive the regressed Y over the whole sample
Yfinal=X(:,[1 C(indexfinal).predictors+1])*Bfinal;

% plot Y_regr vs Y and provide correlation between them
[r p]=corrcoef([Yg Yfinal]);

figure;
plot(Yg,Yfinal,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Yg;Yfinal]) max([Yg;Yfinal])]);
ylim([min([Yg;Yfinal]) max([Yg;Yfinal])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['Optimized regression r= ' num2str(r(1,2),'%2.2f') ' with the predictors ' num2str(C(indexfinal).predictors)]);

