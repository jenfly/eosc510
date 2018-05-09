%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lab 11 (16 Nov 2017) 
% NN modeling of human faces 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This dataset contains grayscale images of 38 individuals and around 64 near frontal images under different illuminations per individual. 
% The faces are proces so that they are cropped and centered, with 32 x 32 pixels each

% NOTE: in order to use SOM methods you need to have
% SOM Toolbox installed from the website:
% http://www.cis.hut.fi/projects/somtoolbox/
% in your Maltab path include the downloaded folder with som toolbox

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

% plot just the first face 
figure;
Y1=reshape(fea(1,:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);
caxis([0 1]);

% select middle column, first 10 pixels (to represent the inout data into NN model)
for j=1:size(fea,1)
matrix=reshape(fea(j,:),[faceH,faceW]);
xdata(:,j)=matrix(1:10,16); 
end

xdata=xdata';

% plot just the first face with NaN
feaNaN(1:size(fea,1),1:size(fea,2))=NaN;
feaNaN(:,15*32+1:15*32+10)=fea(:,15*32+1:15*32+10);

figure;
Y1=reshape(feaNaN(1,:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);
caxis([0 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NN modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for training pick, for example, first first 2000 data points 
% leave the remaining points for forecast (test)
Nin=2000;
Nr=size(fea,1)-Nin;
xdata_in=xdata(1:Nin,:);
xdata_test=xdata(Nin+1:end,:);



