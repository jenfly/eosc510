%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lab 7 with solutions (19 Oct 2017)
% FSA on human face data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This dataset contains grayscale images of 38 individuals and around 64 near frontal images under different illuminations per individual. The faces are processed so that they are cropped and centered, with 32 x 32 pixels each.

load YaleB_32x32.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scale the features (pixel values) to [0,1]
maxValue = max(max(fea));
fea = fea/maxValue;
%===========================================

% plot 100 first images all togeter in one figure
faceW = 32;
faceH = 32;
numPerLine = 10;
ShowLine = 10;

Y = zeros(faceH*ShowLine,faceW*numPerLine);
for i=0:ShowLine-1
   for j=0:numPerLine-1
     Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW) = reshape(fea(i*numPerLine+j+1,:),[faceH,faceW]);
   end
end

figure;
imagesc(Y);
colormap(gray);

% example of a plot (just the first face) 
figure;
subplot(1,2,1)
Y1=reshape(fea(1,:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);

% plot the same figure but in 1-D
x=fea(1,:)'; % convert the first image to vector
N=length(x);  % number of points
L=N;  % lenght of the series
dx=L/N;  % spatial sampling interval

subplot(1,2,2);
plot([1:N],x,'b-');
xlabel('points');
title('image in 1-D');
xlim([1 N]);

% FILTERING WITH A TAPER
% define taper function
taper_function = ones(size(x));
taper_length = 20;
taper_function(1:taper_length) = sin(pi*(0:taper_length-1)/(2*taper_length)).^2;
taper_function(end-taper_length+1:end) = sin(pi*(taper_length-1:-1:0)/(2*taper_length)).^2;

figure;
subplot(2,1,1)
plot([1:N],taper_function,'b-');
xlabel('time');
ylabel('taper function');
xlim([1 N]);

subplot(2,1,2)
plot([1:N],x,'r-',[1:N],taper_function.*x,'b-');
xlabel('time');
ylabel('x2');
legend('original','tapered');
xlim([1 N]);

% tapering
xft_taper = fft(taper_function.*x);

% FFT of x (x(t) --> X(omega)
xft = xft_taper;

% plot the powerspectrum
power_spectrum = abs(xft.^2);
% if we want S to be calculated the same as in periodogram:
% Sm=(2/N*dt)*(x2ft^2); 
S=2*power_spectrum./(N*dx);
S(1)=power_spectrum(1)/(N*dx);
S(1)=0;  % set the mean value to zero

figure; 
subplot(2,1,1)
plot([1:N],S,'b-');
xlabel('points');
ylabel('S of tapered signal');
xlim([1 N]);

% wave number k (same as omega in time domain):
k=[0:N]*2*pi/L;

% plot S versus omega
subplot(2,1,2)
plot(k(1:N/2+1),S(1:N/2+1),'b.-',-k(2:N/2+1),S(end:-1:N/2+1),'r.-'); 
xlabel('angular frequency');
ylabel('S of tapered signal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply low-pass filter (k < 1)
index=170;  % manually identified treshold for k =1
xft_filter = zeros(size(xft));
xft_filter(1:index) = xft(1:index);
xft_filter(end+2-index:end) = xft(end+2-index:end);

power_spectrum = abs(xft_filter.^2);
Sfilter=2*power_spectrum./(N*dx);
Sfilter(1)=0;

figure; 
subplot(1,2,1)
plot([1:N],Sfilter,'b-');
xlabel('points')
ylabel('S of high-pass filter');
xlim([0 N]);

% apply inverse FFT to get back the filtered x(t)
x_filter = ifft(xft_filter);
subplot(1,2,2)
plot([1:N],x,'r-',[1:N],x_filter,'b-');
legend('original','low-pass filtered signal');
xlim([0 N]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply high-pass filter (k > 1)
index=170; % manually identified treshold for k =1
xft_filter = xft;
xft_filter(1:index) = zeros(index,1);
xft_filter(end+2-index:end) = zeros(index-1,1);
xft_filter(1)=xft(1); % keep the power that shows the mean value of the data

power_spectrum = abs(xft_filter.^2);
Sfilter=2*power_spectrum./(N*dx);
Sfilter(1)=0;

figure; 
subplot(1,2,1)
plot([1:N],Sfilter,'b-');
xlabel('points')
ylabel('S of high-pass filter');
xlim([0 N]);

% apply inverse FFT to get back the filtered x(t)
x_filter = ifft(xft_filter);
subplot(1,2,2)
plot([1:N],x,'r-',[1:N],x_filter,'b-');
legend('original','high-pass filtered signal');
xlim([0 N]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q1: apply a filter that keeps only the signals associated with the lowest wavelengths (k small). 
index=20;  % manually identified treshold 
xft_filter = zeros(size(xft));
xft_filter(1:index) = xft(1:index);
xft_filter(end+2-index:end) = xft(end+2-index:end);

power_spectrum = abs(xft_filter.^2);
Sfilter=2*power_spectrum./(N*dx);
Sfilter(1)=0;

figure; 
subplot(1,2,1)
plot([1:N],Sfilter,'b-');
xlabel('points')
ylabel('S of high-pass filter');
xlim([0 N]);

% apply inverse FFT to get back the filtered x(t)
x_filter = ifft(xft_filter);
subplot(1,2,2)
plot([1:N],x,'r-',[1:N],x_filter,'b-');
legend('original','filtered signal');
xlim([0 N]);


