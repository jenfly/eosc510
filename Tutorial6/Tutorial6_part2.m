%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 6 (10 Oct 2017)
% Application of Fourier spectral analysis (FSA)
% on a synthetic data
% this code is the same as Tutorial6_part1.m  but uses 
% the timeseries y that is 310 points long 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% Example 1: autospectrum
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
eng=[0.5 1 2];
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












