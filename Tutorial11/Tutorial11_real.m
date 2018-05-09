%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 11 (14 Nov 2017)
% MLP NN model on real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gridded monthly SST and SLP data for Tropical Pacific
% from ERA Interim reanalysis
% Period: Jan 1979 to Jun 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load SST and SLP 3-month runnng mean anomalies 
load SST_SLP_anomalies.mat;

SSTdata=xdata;
SLPdata=ydata;

clear xdata
clear ydata

months={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};

% start with Feb 1997 and end with May 2015
time=[1979+2/12:1/12:2015+5/12];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA on both datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA on SST:
[eigenvectors_SST,PCs_SST,eigenvalues_SST]=princomp(SSTdata);

% contribution of each mode to total variance:
variance_SST=eigenvalues_SST./sum(eigenvalues_SST);

% plot first 4 modes (eigenvectors and PCs)
figure; 
for i=1:4
subplot(4,2,2*i-1)
var01=eigenvectors_SST(:,i);
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
colorbar;
title(['SST eigenvector var=',num2str(variance_SST(i),'%2.2f')]);

subplot(4,2,2*i)
plot(time,PCs_SST(:,i),'b-')
xlabel('time');
title(['PC',num2str(i)]);
xlim([1979 2015.7]);
end

% PCA on SLP:
[eigenvectors_SLP,PCs_SLP,eigenvalues_SLP]=princomp(SLPdata);

% contribution of each mode to total variance:
variance_SLP=eigenvalues_SLP./sum(eigenvalues_SLP);

% plot first 4 modes (eigenvectors and PCs)
figure; 
for i=1:10
subplot(10,2,2*i-1)
var01=eigenvectors_SLP(:,i);
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
colorbar;
title(['SLP eigenvector var=',num2str(variance_SLP(i),'%2.2f')]);

subplot(10,2,2*i)
plot(time,PCs_SLP(:,i),'b-')
xlabel('time');
title(['PC',num2str(i)]);
xlim([1979 2015.7]);
end

% take first 10 PCs in SLP as x variables in MLP NN model
% to model first PC in SST
% for training pick only first 35 years
% leave the remaining years for forecast (test)
Nin=25*12-1;
xdata_in=PCs_SLP(1:Nin,1:10);
ydata_out=PCs_SST(1:Nin,1);

xdata_test=PCs_SLP(Nin+1:end,1:10);
ydata_test=PCs_SST(Nin+1:end,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model 1 
% MLP NN model using 1 hidden layer with 10 neurons
% and Levenberg-Marquardt optimization 
% run 100 times using random sample for input data (bagging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0  % change to 1 if you want to run this part of the script
no_runs=100;   % note that it takes several minutes to run, so you might test it first with smaller number of runs
number_of_hidden_neurons=10; % here can choose the number of neurons

for kk=1:no_runs
% pick 50 random points from x and ydata
index = randi(Nin,[1 Nin]);
xdata1=xdata_in(index,:);
ydata1=ydata_out(index,:);

% apply MLP NN model
cputime0 = cputime;
net = feedforwardnet(number_of_hidden_neurons,'trainlm');  % here can change to 'trainbr' for BNN
net.trainParam.epochs = 2000;
net.trainParam.show = 500;
[net,tr]= train(net,xdata1',ydata1');
fprintf(1,'\n# cputime = %11.4g\n',cputime-cputime0); cputime0=cputime;

% simulate ouput for xdata_test
ymodel(kk,:) = net(xdata_test');
% find RMSE between this output and sine signal (y)
rmse(kk) = sqrt(mse(ymodel(kk,:)-ydata_test'));
end

% calculate ensamble mean and RMSE
ymodel_mean=mean(ymodel);
rmse(kk+1) = sqrt(mse(ymodel_mean-ydata_test'));

% plot RMSE
figure;
plot([1:no_runs],rmse(1:no_runs),'bo',[no_runs+1],rmse(no_runs+1),'ro');
xlabel('model run');
ylabel('RMSE');

% plot the simulations and ensemble mean
figure;
plot(time(Nin+1:end),ymodel,'y-');
hold on
plot(time(Nin+1:end),ymodel_mean,'b-',time(Nin+1:end),ydata_test,'r-');
xlabel('x');
ylabel('y');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model 2 
% MLP NN model using 1 hidden layer with 1 to 20 neurons
% and Levenberg-Marquardt optimization 
% run 30 times using random sample for input data (bagging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0  % chnage 0 to 1 if you want this part to run
no_runs=30;   % note that it takes several minutes to run
m_max=20;  % maximim number of hidden neurons

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
ymodel_ensemblen(m,:)=net(xdata_in');
rmsen(m) = sqrt(mse(ymodeln(m,:)-ydata_test'));

if kk==1
subplot(4,5,m);
plot(time(Nin+1:end),ymodeln(m,:),'b-',time(Nin+1:end),ydata_test,'r-');
title(['m2=',num2str(m)]);
end

end
% find m with minimum RMSE
[rmse(kk) ind(kk)]=min(rmsen);
ymodel(kk,:)=ymodeln(ind(kk),:);
ymodel_ensemble(kk,:)=ymodel_ensemblen(ind(kk),:);
clear rmsen
end

ymodel_mean=mean(ymodel);
rmse(kk+1) = sqrt(mse(ymodel_mean-ydata_test'));

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

% plot the ensemble and ensemble mean
figure;
plot(time(Nin+1:end),ymodel,'y-');
hold on
plot(time(Nin+1:end),ymodel_mean,'b-',time(Nin+1:end),ydata_test,'r-');
xlabel('x');
ylabel('y');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MLP on the ensemble runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m = 1:m_max  
net = feedforwardnet(m,'trainlm'); % use 'trainbr' for BNN
net.trainParam.epochs = 2000;
net.trainParam.show = 500;
net = train(net,ymodel_ensemble,ydata_out');
ymodeln(m,:) = net(ymodel);
rmse_ensemble(m) = sqrt(mse(ymodeln(m,:)-ydata_test'));
end

% find m with minimum RMSE
[rmse1 ind1]=min(rmse_ensemble);
ymodel_Fensemble=ymodeln(ind1,:);

rmse(kk+2) = sqrt(mse(ymodel_Fensemble-ydata_test'));

figure;
plot([1:no_runs],rmse(1:no_runs),'bo',[no_runs+1],rmse(no_runs+1),'ro',[no_runs+2],rmse(no_runs+2),'go');
xlabel('model run');
ylabel('RMSE');
xlim([0 no_runs+2])

% plot the ensemble, ensemble mean, MLP ensemble mean
figure;
plot(time(Nin+1:end),ymodel,'y-');
hold on
plot(time(Nin+1:end),ymodel_mean,'b-',time(Nin+1:end),ymodel_Fensemble,'g-',time(Nin+1:end),ydata_test,'r-');
xlabel('x');
ylabel('y');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model 3
% Multiple-linear regression model (stepwise) just for comparison 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0  % chnage to 1 if you want this to run

% MLR with original variables
n=size(ydata_out,1);
n2=size(ydata_test,1);
X=[ones(n,1) xdata_in];   % first column needs to consist of ones
Xtest=[ones(n2,1) xdata_test];

b1=regress(ydata_out,X)
y_regr1=X*b1;
y_model1=Xtest*b1;

% stepwise regression with original values
[a SE PVAL INMODEL STATS NEXTSTEP HISTORY]=stepwisefit(xdata_in,ydata_out,'penter',0.01);
% regression coefficients
a=a
% constant coefficient 
a0=STATS.intercept

% take out only the relevant predictors (from stepwisefit)
[ind1 ind2]=find(INMODEL == 1);
b=regress(ydata_out,X(:,[1 ind2+1]))
y_regr=X(:,[1 ind2+1])*b;
y_model=Xtest(:,[1 ind2+1])*b;

SSR=sum((y_model-mean(ydata_test)).^2);
SST=sum((ydata_test-mean(ydata_test)).^2);
R2=SSR/SST;

xline=[min([ydata_test;y_model]):max([ydata_test;y_model])]; 
figure;
plot(ydata_test,y_model,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([ydata_test;y_model]) max([ydata_test;y_model])]);
ylim([min([ydata_test;y_model]) max([ydata_test;y_model])]);
xlabel('y');
ylabel('y_{regr}')
title(['Stepwise regression R^2= ' num2str(R2,'%2.2f')]);

rmse1= sqrt(mse(y_model-ydata_test));

% plot the simulations 
figure;
subplot(2,1,1)
plot(time(Nin+1:end),y_model,'b-',time(Nin+1:end),ydata_test,'r-');
label('MLR','original','MLP NN')
xlabel('time');
ylabel('SST PC1');
title('stepwise');

subplot(2,1,2)
plot(time(Nin+1:end),y_model1,'b-',time(Nin+1:end),ydata_test,'r-');
label('time');
ylabel('SST PC1');
title('MLR');

end

