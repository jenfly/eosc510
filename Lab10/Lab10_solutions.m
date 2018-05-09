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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NN modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% take first 5 PCs in x variables in MLP NN model
% for training pick first 90 points
% leave the remaining points for forecast (test)

Nin=90
xdata_in=PCs(1:Nin,1:5);
ydata_out=Yg(1:Nin,1);

xdata_test=PCs(Nin+1:end,1:5);
ydata_test=Yg(Nin+1:end,1);

xdata=PCs(:,1:5);
ydata=Yg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run 1 
% MLP NN model using 1 hidden layer with 10 neuronsc
% and Levenberg-Marquardt optimization 
% run N times using random sample for input data (bagging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1  % change to 1 if you want to run this part of the script

no_runs=50;   % note that it takes several minutes to run, so you might test it first with smaller number of runs
number_of_hidden_neurons=10; % here can choose the number of neurons

for kk=1:no_runs
% pick 50 random points from x and ydata
index = randi(Nin,[1 Nin]);
xdata1=xdata_in(index,:);
ydata1=ydata_out(index,:);

% apply MLP NN model
cputime0 = cputime;
net = feedforwardnet(number_of_hidden_neurons,'trainlm');  
net.trainParam.epochs = 2000;
net.trainParam.show = 500;
[net,tr]= train(net,xdata1',ydata1');
fprintf(1,'\n# cputime = %11.4g\n',cputime-cputime0); cputime0=cputime;

% simulate ouput for xdata_test
ymodel(kk,:) = net(xdata_test');
% find RMSE between this output and test data
rmse(kk) = sqrt(mse(ymodel(kk,:)-ydata_test'));
end

% calculate ensamble mean and RMSE
ymodel_mean=mean(ymodel);
rmse(kk+1) = sqrt(mse(ymodel_mean-ydata_test'));

RMSE(1)=rmse(kk+1);

% plot RMSE
figure;
plot([1:no_runs],rmse(1:no_runs),'bo',[no_runs+1],rmse(no_runs+1),'ro');
xlabel('model run');
ylabel('RMSE');

[r p]=corrcoef([ydata_test ymodel' ymodel_mean']);

figure;
plot(ydata_test,ymodel_mean,'bo'); hold on
plot(ydata_test,ydata_test,'k-','LineWidth',1);
xlabel('g obs');
ylabel('g model')
title(['MLP model mean w/ 10 hidden neurons r= ' num2str(r(1,end),'%2.2f')]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run 2 
% MLP NN model using 1 hidden layer with 1 to 10 neurons
% and Levenberg-Marquardt optimization 
% run N times using random sample for input data (bagging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1  % change 0 to 1 if you want this part to run

no_runs=30;   % note that it takes several minutes to run
m_max=10;  % maximim number of hidden neurons

for kk=1:no_runs

kk=kk
index = randi(Nin,[1 Nin]);
xdata1=xdata_in(index,:);
ydata1=ydata_out(index,:);

for m = 1:m_max  
net = feedforwardnet(m,'trainlm'); % use 'trainbr' for BNN
net.trainParam.epochs = 2000;
net.trainParam.show = 500;
net = train(net,xdata1',ydata1');
ymodeln(m,:) = net(xdata_test');
rmsen(m) = sqrt(mse(ymodeln(m,:)-ydata_test'));

end

% find m with minimum RMSE
[rmse(kk) ind(kk)]=min(rmsen);
ymodel(kk,:)=ymodeln(ind(kk),:);
clear rmsen
end

ymodel_mean=mean(ymodel);
rmse(kk+1) = sqrt(mse(ymodel_mean-ydata_test'));
RMSE(2)=rmse(kk+1);

% plot RMSE
figure;
subplot(2,1,1);
plot([1:no_runs],rmse(1:no_runs),'bo',[no_runs+1],rmse(no_runs+1),'ro');
xlabel('model run');
ylabel('RMSE');
xlim([0 no_runs+1])

subplot(2,1,2);
plot([1:no_runs],ind,'ko');
xlabel('model run');
ylabel('number of hidden neurons');
xlim([0 no_runs+1])


[r p]=corrcoef([ydata_test ymodel_mean']);

figure;
plot(ydata_test,ymodel_mean,'bo'); hold on
plot(ydata_test,ydata_test,'k-','LineWidth',1);
xlabel('g obs');
ylabel('g model')
title(['MLP model mean w/ optimal hidden neurons r= ' num2str(r(1,end),'%2.2f')]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional: MLR and stepwise on first five modes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MLP with standardized variables
X=[ones(Nin,1) xdata_in];   
Xm=[ones(length(ydata_test),1) xdata_test];
[b2,BINT2,R2,RINT2,STATS2]=regress(ydata_out,X);
y_regrMLR=Xm*b2;

RMSE(3) = sqrt(mse(y_regrMLR-ydata_test));

% stepwise regression with standardized values
[a2 SE PVAL INMODEL STATS NEXTSTEP HISTORY]=stepwisefit(X(:,2:end),ydata_out,'penter',0.05);

% regression coefficients
a2=a2

% constant coefficient 
a20=STATS.intercept

% take out only the relevant predictors (from stepwisefit)
[ind1 ind2]=find(INMODEL == 1);
b=regress(ydata_out,X(:,[1 ind2+1]))
y_regrSTEP=Xm(:,[1 ind2+1])*b;

RMSE(4) = sqrt(mse(y_regrSTEP-ydata_test));

[r p]=corrcoef([ydata_test y_regrMLR y_regrSTEP]);

figure;
subplot(2,1,1)
plot(ydata_test,y_regrMLR,'bo'); hold on
plot(ydata_test,ydata_test,'k-','LineWidth',1);
xlabel('Y');
ylabel('Y_{regr}')
title(['MLR r= ' num2str(r(1,2),'%2.2f')]);

subplot(2,1,2)
plot(ydata_test,y_regrSTEP,'bo'); hold on
plot(ydata_test,ydata_test,'k-','LineWidth',1);
xlabel('Y');
ylabel('Y_{regr}')
title(['Stepwise regression r= ' num2str(r(1,3),'%2.2f')]);
