%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lab 8 (24 Oct 2017) 
% Clustering of images in PC space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This dataset contains grayscale images of 38 individuals and around 64 near frontal images under different illuminations per individual. 

% The faces are proces so that they are cropped and centered, with 32 x 32 pixels each

%The code is from Lab3:
% performs PCA on the data (finds eigenvectors, PCs) 

load YaleB_32x32.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scale the features (pixel values) to [0,1]
maxValue = max(max(fea));
fea = fea/maxValue;
%===========================================

% plot 100 first images all togeter on one figure
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
title('first 100 faces');

% plot just the first face 
figure;
Y1=reshape(fea(1,:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);
title('first face');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[eigenvectors,PCs,eigenvalues]=princomp(fea);




