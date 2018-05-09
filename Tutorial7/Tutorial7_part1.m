%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 7 (17 Oct 2017)
% Filtering on synthetic data
% this code deals with the data that are 300 points long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
% Example 1: filtering in the frequency domain
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

N=300;  % number of points
T=300;  % length of the record
dt=T/N; % time interval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: FFT of x2 (x(t) --> X(omega)
x2ft = fft(x2);

% plot the powerspectrum
power_spectrum = abs(x2ft.^2);
% if we want S to be calculated the same as in periodogram:
% Sm=(2/N*dt)*(x2ft^2); 
S=2*power_spectrum./(N*dt);
S(1)=power_spectrum(1)/(N*dt);

figure; 
plot(S);
xlabel('points');
ylabel('S');

omega=[0:N]*2*pi/T;
% plot S versus omega
figure; 
plot(omega(1:T/2+1),S(1:T/2+1),'b.-',-omega(2:T/2+1),S(end:-1:T/2+1),'r.-'); 
xlabel('angular frequency');
ylabel('S');
xlim([-1 1]);

% Step 2: apply filter on X(omega)to retrieve the signal with lowest frequency 
% the index for the power peak at the lowest fequency is 3
index=3;
x2ft_filter = zeros(size(x2ft));
x2ft_filter(index) = x2ft(index);
x2ft_filter(end+2-index) = x2ft(end+2-index);

figure; 
plot(abs(x2ft_filter).^2);
xlabel('points')
ylabel('filtered spectrum');

% Step 3: apply inverse FFT to get back the filtered x(t)
x2_filter = ifft(x2ft_filter);
figure;
subplot(3,1,1)
plot(t,y(1,:),'r-',t,x2_filter,'b-');
legend('original','FSA');

% retrieve the signal with medium frequency
index=5;
x2ft_filter = zeros(size(x2ft));
x2ft_filter(index) = x2ft(index);
x2ft_filter(end+2-index) = x2ft(end+2-index);

x2_filter = ifft(x2ft_filter);
subplot(3,1,2)
plot(t,y(2,:),'r-',t,x2_filter,'b-');
legend('original','FSA');

% retrieve the signal with highest frequency
index=16;
x2ft_filter = zeros(size(x2ft));
x2ft_filter(index) = x2ft(index);
x2ft_filter(end+2-index) = x2ft(end+2-index);

x2_filter = ifft(x2ft_filter);
subplot(3,1,3)
plot(t,y(3,:),'r-',t,x2_filter,'b-');
legend('original','FSA');




