%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lab 6 with solutions (12 Oct 2017)
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
Y1=reshape(fea(1,:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);

% plot the same figure but in 1-D
figure;
plot(fea(1,:));
title('first face in 1-D');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q1: What is the spectrum of this 1-D face? 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=fea(1,:)'; % convert the first image to vector
N=length(x);  % number of points
dx=1;  % spatial sampling interval
L=N;  % lenght of the series

% option 1:
% Note: wave number k is same as omega in time domain:
[Sx1, k1] = periodogram(fea(1,:)',[],N,1);

% option 2:
xfft = fft(x);
Sx2n = abs(xfft.^2);

% if we want S to be calculated the same as in periodogram function then:
Sx2=2*Sx2n./(N*dx);
Sx2(1)=Sx2n(1)./(N*dx);

% wave number k (same as omega in time domain):
k2=[0:N]*2*pi/L;

figure; 
plot(2*pi*k1,Sx1,'b-',k2(1:N/2+1),Sx2(1:N/2+1),'-');
xlabel('k');
ylabel('spectrum of 1-D face');
%xlim([0 2]);

% there is huge energy on k=0 (= mean value of the data); 
% lets ignore the mean value by setting Sx(1) to zero, and only look the spectrum for k>0 
Sx2(1)=0;
Nx=32;
Ny=32;
dx=1;
dy=1;
Lx=32;
Ly=32;

figure; 
plot(k2(1:N/2+1),Sx2(1:N/2+1),'-');
xlabel('k');
ylabel('spectrum of 1-D face');
%xlim([0 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q2: What is 2-D spectrum of the 2-D face (first face in the dataset)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now lets do fft in 2D:
Y1=reshape(fea(1,:),[faceH,faceW]);
Sxy = fft2(Y1);
Sxy = 4*abs(Sxy.^2)./(Nx*Ny*dx*dy);

% lets ignore the mean value by setting Sxy(1,1) to zero, and only look the spectrum for kx>0 ky>0
Sxy(1,1) = 0;

% lets center the spectrum so that Sxy(1,1) is at the center of the image
Sxy = fftshift(fftshift(Sxy,1),2);

kx=[-Nx/2:Nx/2-1]*2*pi/Lx;
ky=[-Ny/2:Ny/2-1]*2*pi/Ly;

figure; 
imagesc(kx,ky,Sxy) 
colormap jet
colorbar
title('2-D spectrum of the 1st face');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q3: What is 2-D spectrum of the 2-D face plotted below (77th face in the dataset)? 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrum of the 77th face
Y1=reshape(fea(77,:),[faceH,faceW]);
Sxy = fft2(Y1);
Sxy = 4*abs(Sxy.^2)./(Nx*Ny*dx*dy);

% lets ignore the mean value by setting Sxy(1,1) to zero, and only look the spectrum for kx>0 ky>0
Sxy(1,1) = 0;

% lets center the spectrum so that Sxy(1,1) is at the center of the image
Sxy = fftshift(fftshift(Sxy,1),2);

figure; 
imagesc(kx,ky,Sxy) 
colormap jet
colorbar
title('2-D spectrum of the 77th face');

% we can see that here the high wevenumber (shorth wavelengths) dominate the power spectrum -> this signal is
% noise and not the actual shading of the face

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q4 (optional): Calculate the 2-D spectrum of all the faces in the dataset. 
% What are the most characteristic spectra in the dataset? 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:size(fea,1)
Y1=reshape(fea(j,:),[faceH,faceW]);
Sxy = fft2(Y1);
Sxy = 4*abs(Sxy.^2)./(Nx*Ny*dx*dy);
Sxy(1,1) = 0;
Sxy = fftshift(fftshift(Sxy,1),2);
% convert 2-D Sxy to 1-D column and save it under spectra for that column:
spectra(j,:)=reshape(Sxy,1,size(fea,2));
end

% performing PA on spectra:
[eigenvectors,PCs,eigenvalues]=princomp(spectra);

% contribution of each mode to total variance:
variance=eigenvalues./sum(eigenvalues);

% plot the variance for each mode
figure; 
subplot(2,1,1)
plot([1:length(variance)],variance,'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

subplot(2,1,2)
plot([1:10],variance(1:10),'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

% we can see that the first mode carries almost all the variance

% plot the first 3 modes (eigenvectors and PCs)
% but actually we care only about the first mode
figure; 
modes=3;
for j=1:modes
subplot(2,modes,j)
Y1=reshape(eigenvectors(:,j),[faceH,faceW]);
imagesc(kx,ky,Y1);
colormap(jet);
title(['e_',num2str(j)]);

subplot(2,modes,j+modes)
plot(PCs(:,j));
xlabel('index');
title(['PC',num2str(j)]);
xlim([1 size(fea,1)])
end







