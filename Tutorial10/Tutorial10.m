%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 10 (7 Nov 2017)
% MLP NN modeling on a synthetic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating the input data
x=[1:2:100];
k=2*pi/50;
y=sin(k*x);

figure;
plot(x,y,'b-');
title('sine signal');

% generate training data w/ noise
load noise.mat;
noise1=noise(1,[1:2:100]);
factor=std(y)/std(noise1);

ydata=sin(k*x)+factor*noise1;

figure;
plot(x,ydata,'bo',x,y,'r-');
title('sine signal with noise');
legend('training data','sine signal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1 
% MLP NN model using 1 hidden layer with 10 neuronsc
% and Levenberg-Marquardt optimization 
% run N times using random sample for input data (bagging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1  % chnage to 0 if you don't want this example to run

N=30; % number of model runs
for kk=1:N
% pick 50 random points from x and ydata
index = randi(50,[1 50]);
xdata1=x(index);
ydata1=ydata(index);

% apply MLP NN model
cputime0 = cputime;
net = feedforwardnet(10,'trainlm');
%net=fitnet(10,'trainlm'); % this is equivalent command 
net.trainParam.epochs = 2000;
net.trainParam.show = 500;
net= train(net,xdata1,ydata1);
fprintf(1,'\n# cputime = %11.4g\n',cputime-cputime0); cputime0=cputime;

% simulate ouput for x
ymodel(kk,:) = net(x);
% find RMSE between this output and sine signal (y)
rmse(kk) = sqrt(mse(ymodel(kk,:)-y));
end

% calculate ensamble mean and RMSE
ymodel_mean=mean(ymodel);
rmse(kk+1) = sqrt(mse(ymodel_mean-y));

% plot RMSE
figure;
plot([1:N],rmse(1:N),'bo',[N+1],rmse(N+1),'ro');
xlabel('model run');
ylabel('RMSE');

% plot the simulations and ensemble mean
figure;
plot(x,ydata,'ko',x,ymodel,'k-');
hold on
plot(x,ymodel_mean,'b-',x,y,'r-');
xlabel('x');
ylabel('y');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2 
% MLP NN model using 1 hidden layer with 1 to 12 neurons
% and Levenberg-Marquardt optimization 
% run N times using random sample for input data (bagging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0  % change 0 to 1 if you want this part to run

m_max=12;
N=20; % number of model runs

for kk=1:N
kk=kk
index = randi(50,[1 50]);
xdata1=x(index);
ydata1=ydata(index);

for m = 1:m_max  
net = feedforwardnet(m,'trainlm');
net.trainParam.epochs = 2000;
net.trainParam.show = 500;
net = train(net,xdata1,ydata1);
ymodeln(m,:) = net(x);
rmsen(m) = sqrt(mse(ymodeln(m,:)-y));

if kk==1
subplot(3,4,m);
plot(x,ymodeln(m,:),'b-');
title(['m2=',num2str(m)]);
end

end
% find m with minimum RMSE
[rmse(kk) ind(kk)]=min(rmsen);
ymodel(kk,:)=ymodeln(ind(kk),:);
clear rmsen
end

ymodel_mean=mean(ymodel);
rmse(kk+1) = sqrt(mse(ymodel_mean-y));

% plot RMSE
figure;
subplot(2,1,1);
plot([1:N],rmse(1:N),'bo',[N+1],rmse(N+1),'ro');
xlabel('model run');
ylabel('RMSE');
xlim([0 N+2]);

subplot(2,1,2);
plot([1:N],ind,'ko');
xlabel('model run');
ylabel('number of hidden neurons');
xlim([0 N+2]);

% plot the ensemble and ensemble mean
figure;
plot(x,ydata,'ro',x,ymodel,'k-');
hold on
plot(x,ymodel_mean,'b-',x,y,'r-');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 3 
% MLP BNN model using 1 hidden layer with 1 to 12 neurons
% run 10 times using random sample for input data (bagging)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0 % change 0 to 1 if you want this part to run
% NOTE: this one runs for long

N=10;   % number of model runs 

m_max=12;

for kk=1:N
kk=kk
index = randi(50,[1 50]);
xdata1=x(index);
ydata1=ydata(index);

for m = 1:m_max  
net = feedforwardnet(m,'trainbr');
net.trainParam.epochs = 2000;
net.trainParam.show = 500;
net = train(net,xdata1,ydata1);
ymodeln(m,:) = net(x);
rmsen(m) = sqrt(mse(ymodeln(m,:)-y));

if kk==1
subplot(3,4,m);
plot(x,ymodeln(m,:),'b-');
title(['m2=',num2str(m)]);
end

end
% find m with minimum RMSE
[rmse(kk) ind(kk)]=min(rmsen);
ymodel(kk,:)=ymodeln(ind(kk),:);
clear rmsen
end

ymodel_mean=mean(ymodel);
rmse(kk+1) = sqrt(mse(ymodel_mean-y));

% plot RMSE
figure;
subplot(2,1,1);
plot([1:N],rmse(1:N),'bo',[N+1],rmse(N+1),'ro');
xlabel('model run');
ylabel('RMSE');
xlim([0 N+2]);

subplot(2,1,2);
plot([1:N],ind,'ko');
xlabel('model run');
ylabel('number of hidden neurons');
xlim([0 N+2]);

% plot the ensemble and ensemble mean
figure;
plot(x,ydata,'ro',x,ymodel,'k-');
hold on
plot(x,ymodel_mean,'b-',x,y,'r-');
end
