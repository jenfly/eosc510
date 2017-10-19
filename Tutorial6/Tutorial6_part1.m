%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 6 (10 Oct 2017)
% Application of Fourier spectral analysis (FSA)
% on a synthetic data
% this code deals with the data that are 300 points long (from  Tutorial)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% Example 1: autospectrum
%%%%%%%%%%%%%%%%%% 
% creating x1 and x2 as a superposition of synthetic signals 
k=2*pi/50;
omega=[2*pi/150 2*pi/75 2*pi/20]; % = [0.0419    0.0838    0.3142]
x=[1:100];
t=[1:300];

y(1,:)=0.5*sin(k*x(1)-omega(1)*t);
y(2,:)=sin(k*x(17)-omega(2)*t);
y(3,:)=2*sin(k*x(70)-omega(3)*t);
y(4,:)=0.01*t-1;

x1=sum(y);

% subtract linear regression line
tn(1:300)=1;
tn=[tn' t'];
b=regress(x1',tn);
trend=tn*b;

x2=x1-trend';	

figure;
subplot(3,1,1)
plot(y');
xlabel('time');

subplot(3,1,2)
plot(x1);
xlabel('time');
ylabel('x1');

subplot(3,1,3)
plot(x2);
xlabel('time');
ylabel('x2');

% calulating autospectrum S
[Px1, w1] = periodogram(x1,[],size(x1,2));
[Px2, w2] = periodogram(x2,[],size(x2,2));

% plotting autospectrum 
figure;
subplot(2,1,1)
plot(w1,Px1);
title('Spectrum of x1');
xlabel('angular frequency');
ylabel('spectrum');

subplot(2,1,2)
plot(w2,Px2);
title('Spectrum of x2');
xlabel('angular frequency');
ylabel('spectrum');

% calulating autospectrum S using different option in matlab
% the one with default (2^n) number of point used in the FFT
[Px1n, w1n] = periodogram(x1);
[Px2n, w2n] = periodogram(x2);

% plotting autospectrum 
figure;
subplot(2,1,1)
plot(w1n,Px1n);
title('Spectrum of x1');
xlabel('angular frequency');
ylabel('spectrum');

subplot(2,1,2)
plot(w2n,Px2n);
title('Spectrum of x2');
xlabel('angular frequency');
ylabel('spectrum');

% variance in Fourier bands (% of total)
[energy1 ind1] = sort(100*Px1/sum(Px1),'descend');
[energy2 ind2] = sort(100*Px2/sum(Px2),'descend');

% angular frequencies for first 10 peaks
freq1=w1(ind1(1:10));
freq1n=w1n(ind1(1:10));

freq2=w2(ind2(1:10));
freq2n=w2n(ind2(1:10));

% real energy(%) in the data
% amplitude^2
eng=[0.5 1 2];  % these are amplitudes of the known signals (y1, y2, y3)
eng=100*eng.^2./sum(eng.^2);

% plot energy per frequency (real and derived from FSA)
figure;
subplot(2,1,1)
plot(omega,eng,'bo',freq1,energy1(1:10),'r*');
xlabel('omega');
ylabel('x1 energy');

subplot(2,1,2)
plot(omega,eng,'bo',freq2,energy2(1:10),'r*');
xlabel('omega');
ylabel('x2 energy');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for aliasing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lets assume we sample the signal x1 at fewer points
% the signal with max frequancy has a time period (theta) = 20
% we need at least 2 points to cover this period 
% (the original signal has 20 point in theta, do dt=1)

% set dt;
dN=3; % new sampling interval 
ind=[1:dN:300];
dt=300/length(ind);
x2new=x2(ind);
freq=1/dt; % new sampling frequency
omega_N=pi/dt; % Nyquist angular frequency

figure;
subplot(2,1,1)
plot(t,x2,'b.-'); hold on
plot(ind,x2new,'ro');
xlabel('time');
ylabel('x2');

% calulating autospectrum S
[Px2, w2] = periodogram(x2,[],length(t),1);
[Px2new, w2new] = periodogram(x2new,[],length(ind),freq);

om(1:2)=omega_N;
yy=[0 max(Px2new)];
% plotting autospectrum 
subplot(2,1,2)
plot(2*pi*w2,Px2,'bo-',2*pi*w2new,Px2new,'r*-'); hold on
plot(om,yy,'k--','LineWidth',2);
xlabel('angular frequency');
ylabel('spectrum');
legend('original','new','Nyquist freq.');
xlim([0 1.1]);

% set dt;
dN=9; % new sampling interval 
ind=[1:dN:300];
dt=300/length(ind);
x2new=x2(ind);
freq=1/dt; % new sampling frequency
omega_N=pi/dt; % Nyquist angular frequency

figure;
subplot(2,1,1)
plot(t,x2,'b.-'); hold on
plot(ind,x2new,'ro');
xlabel('time');
ylabel('x2');

% calulating autospectrum S
[Px2, w2] = periodogram(x2,[],length(t),1);
[Px2new, w2new] = periodogram(x2new,[],length(ind),freq);

om(1:2)=omega_N;
yy=[0 max(Px2new)];
% plotting autospectrum 
subplot(2,1,2)
plot(2*pi*w2,Px2,'bo-',2*pi*w2new,Px2new,'r*-'); hold on
plot(om,yy,'k--','LineWidth',2);
xlabel('angular frequency');
ylabel('spectrum');
legend('original','new','Nyquist freq.');
xlim([0 0.5]);

% set dt;
dN=12; % new sampling interval 
ind=[1:dN:300];
dt=300/length(ind);
x2new=x2(ind);
freq=1/dt; % new sampling frequency
omega_N=pi/dt; % Nyquist angular frequency

figure;
subplot(2,1,1)
plot(t,x2,'b.-'); hold on
plot(ind,x2new,'ro');
xlabel('time');
ylabel('x2');

% calulating autospectrum S
[Px2, w2] = periodogram(x2,[],length(t),1);
[Px2new, w2new] = periodogram(x2new,[],length(ind),freq);

om(1:2)=omega_N;
yy=[0 max(Px2new)];
% plotting autospectrum 
subplot(2,1,2)
plot(2*pi*w2,Px2,'bo-',2*pi*w2new,Px2new,'r*-'); hold on
plot(om,yy,'k--','LineWidth',2);
xlabel('angular frequency');
ylabel('spectrum');
legend('original','new','Nyquist freq.');
xlim([0 0.5]);








