%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 7 (17 Oct 2016)
% Application of moving-average and singular spectrum analysis (SSA)
% on synthetic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% Moving average 
%%%%%%%%%%%%%%%%%% 
% creating x1 and x2 as a superposition of synthetic signals 
k=2*pi/50;
omega=[2*pi/150 2*pi/75 2*pi/20];
x=[1:100];
t=[1:310];

y(1,:)=0.5*sin(k*x(1)-omega(1)*t);
y(2,:)=sin(k*x(17)-omega(2)*t);
y(3,:)=2*sin(k*x(70)-omega(3)*t);
y(4,:)=0.01*t-1;

x1=sum(y);

% subtract linear regression line
tn(1:310)=1;
tn=[tn' t'];
b=regress(x1',tn);
trend=tn*b;

x2=x1-trend';	

figure;
subplot(3,1,1)
plot(y');
xlabel('time');
title('synthetic signals');

subplot(3,1,2)
plot(x1);
xlabel('time');
ylabel('x1');

subplot(3,1,3)
plot(x2);
xlabel('time');
ylabel('x2');

% define taper function for FSA
taper_function = ones(size(x2));
taper_length = 20;
taper_function(1:taper_length) = sin(pi*(0:taper_length-1)/(2*taper_length)).^2;
taper_function(end-taper_length+1:end) = sin(pi*(taper_length-1:-1:0)/(2*taper_length)).^2;

% Apply running mean (moving average) as a filter in time domain
% from autospectrum we read the frequncies for the peaks
[P, w] = periodogram(taper_function.*x2,[],size(x2,2));
omeg=w([3 5 16 17]);

% calculate window sizes for the moving average from these 
% frequencies from autospectrum
window1=2*pi/omeg(1);
window2=2*pi/omeg(2);
window3=2*pi/mean(omeg(3:4));

% start by applying the smallest window (window 3)
% i.e. removing high-frequency signal

% 1) apply window3 running mean on x2
window=floor(window3);
for i=1:length(x2)-window+1;
    runmean(i)=nanmean(x2(i:i+window-1));
end

x2_runmean=runmean;

% original series - running mean = signal with higest freqency
y3_rec=x2(floor(window3/2):end-floor(window3/2))-x2_runmean;
figure; 
subplot(3,1,1)
title('results from the application of moving-average');
plot(t,y(3,:),'r-',t(floor(window3/2):end-floor(window3/2)),y3_rec,'b-');
legend('original','filtered')

% 2) apply window2 running mean on x2_runmean
clear runmean;

window=floor(window2);
for i=1:length(x2_runmean)-window+1;
    runmean(i)=nanmean(x2_runmean(i:i+window-1));
end

x3_runmean=runmean;

% original series - running mean = signal with higest freqency
t_new=t(floor(window3/2):end-floor(window3/2));
y2_rec=x2_runmean(floor(window2/2):end-floor(window2/2)-1)-x3_runmean;
subplot(3,1,2)
plot(t,y(2,:),'r-',t_new(floor(window2/2):end-floor(window2/2)-1),y2_rec,'b-');

% x3_runmean should resemble y(1,:)
subplot(3,1,3)
plot(t,y(1,:),'r-',t_new(floor(window2/2):end-floor(window2/2)-1),x3_runmean,'b-');

figure;
subplot(3,1,1)
plot(t,x2,'k-',t(floor(window3/2):end-floor(window3/2)),x2_runmean,'g-');
xlabel('x2');
legend('timeseries','running mean');
title('steps in the application of moving average');

subplot(3,1,2)
plot(t(floor(window3/2):end-floor(window3/2)),x2_runmean,'k-',t_new(floor(window2/2):end-floor(window2/2)-1),x3_runmean,'g-');
xlabel('x2 step1');

subplot(3,1,3)
plot(t_new(floor(window2/2):end-floor(window2/2)-1),x3_runmean,'k-');
xlabel('x2 step2');


%%%%%%%%%%%%%%%%%%
% SSA
%%%%%%%%%%%%%%%%%% 
% creating lagged copies of the time series
L=151;
n=size(x1,2);

x1_lag(1:L,1:n-L+1)=NaN;  % reserving space
x2_lag=x1_lag;
for i=1:L
x1_lag(i,:)=x1(1,i:n-L+i);
x2_lag(i,:)=x2(1,i:n-L+i);
end

% PCA applied on x1_lag (w/ trend) and x2_lag (detrended)
[eigenvectors1,PCs1,eigenvalues1]=princomp(x1_lag');
[eigenvectors2,PCs2,eigenvalues2]=princomp(x2_lag');

% contribution of each mode pair to total variance 
variance1=eigenvalues1([1:2:L-1])+eigenvalues1([2:2:L-1]);
variance1=100*variance1./sum(eigenvalues1); % in percent

variance2=eigenvalues2([1:2:L-1])+eigenvalues2([2:2:L-1]);
variance2=100*variance2./sum(eigenvalues2); % in percent

% plotting first 4 pairs of eigenvectors and PCs for x1 (w/ trend)
figure;
for i=1:4
subplot(4,2,i*2-1)
plot(eigenvectors1(:,i*2-1:i*2))
title('SSA eigenvectors for x1');
title(['x1 eigenvectors ',num2str(2*i-1),' ',num2str(2*i)])

subplot(4,2,i*2)
plot(PCs1(:,i*2-1:i*2));
title('PCs for x1');
end

% plotting first 4 pairs of eigenvectors and PCs for x2 (detrended)
figure;
for i=1:4
subplot(4,2,i*2-1)
plot(eigenvectors2(:,i*2-1:i*2))
title('SSA eigenvectors for x2');
title(['x2 eigenvectors ',num2str(2*i-1),' ',num2str(2*i)])

subplot(4,2,i*2)
plot(PCs2(:,i*2-1:i*2));
title('PCs for x1');
end

% reconstructing the three signals from the pair of modes
% The first signal the one with higest frequencies, so this on is reconstructred from 
% the first pair of modes (mode 1 and mode 2)
% The second signal (medium feequancy) is reconstructed from modes 3 and 4
% the third signal (low frequncy) is reconstructed from modes 5 and 6
for i=1:3
x1_lag_rec(i,:)=PCs1(1,i*2-1).*eigenvectors1(:,i*2-1)+PCs1(1,i*2).*eigenvectors1(:,i*2);
x2_lag_rec(i,:)=PCs2(1,i*2-1).*eigenvectors2(:,i*2-1)+PCs2(1,i*2).*eigenvectors2(:,i*2);
end

figure;
subplot(3,1,1)
plot(t,y(3,:),'r-',t(1:151),x1_lag_rec(1,:),'b-',t(1:151),x2_lag_rec(1,:),'g-');
xlabel('time');
xlim([1 151])
ylim([-2.2 2.2]);
legend('original','x1 reconstr','x2 reconstr');
title('reconstruction of signals with SSA');

subplot(3,1,2)
plot(t,y(2,:),'r-',t(1:151),x1_lag_rec(2,:),'b-',t(1:151),x2_lag_rec(2,:),'g-');
xlabel('time');
xlim([1 151])
ylim([-2.2 2.2]);

subplot(3,1,3)
plot(t,y(1,:),'r-',t(1:151),x1_lag_rec(3,:),'b-',t(1:151),x2_lag_rec(3,:),'g-');
xlabel('time');
xlim([1 151])
ylim([-2.2 2.2]);




