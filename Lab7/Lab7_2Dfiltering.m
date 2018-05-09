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

% FILTERING
Nx=32;
Ny=32;
dx=1;
dy=1;
Lx=32;
Ly=32;

figure;
Y1=reshape(fea(77,:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);

% Spectrum of the 77th face
Y1=reshape(fea(77,:),[faceH,faceW]);
Sxy = fft2(Y1);
xy_fft=Sxy;

Sxy = 4*abs(Sxy.^2)./(Nx*Ny*dx*dy);

% lets ignore the mean value by setting Sxy(1,1) to zero, and only look the spectrum for kx>0 ky>0
Sxy(1,1) = 0;

% lets center the spectrum so that Sxy(1,1) is at the center of the image
Sxy_shift = fftshift(fftshift(Sxy,1),2);

kx=[-Nx/2:Nx/2-1]*2*pi/Lx;
ky=[-Ny/2:Ny/2-1]*2*pi/Ly;

figure; 
imagesc(kx,ky,Sxy_shift) 
colormap jet
colorbar
title('shifted 2-D spectrum of the 77th face');

% the actual spectrum (not shifted)
figure; 
imagesc(Sxy) 
colormap jet
colorbar
title('2-D spectrum of the 77th face');

% manually selected points where the noise is in the spectrum
xy_filter=Sxy;
xy_filter([20:27],[4:8])=0;
xy_filter([7:14],[26:30])=0;

figure;
imagesc(xy_filter);
title('image of the filtered spectrum');

xyft_filter=xy_fft;
xyft_filter([20:27],[4:8])=0;
xyft_filter([7:14],[26:30])=0;


% apply inverse FFT to get back the filtered x(t)
xy_new = ifft2(xyft_filter);

figure;
imagesc(xy_new);
colormap(gray);




