%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lab 9 (2 Nov 2017) 
% SOMs of human faces 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOM algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=double(fea);

% chose size of the map for SOM
% I put this in a loop so that I save the quantization error and topographic error 
% for each SOM size 
% SOM sizes I choose: 3 x 3, 4 x 3, 4 x 4, 5 x 3, 5 x 4, 5 x 5
ny_arrey=[3 4 5 4 5 5];
nx_arrey=[3 3 3 4 4 5];

for k=1:length(ny_arrey)

ny_som=ny_arrey(k); nx_som=nx_arrey(k);

en=ny_som*nx_som;
msize=[ny_som nx_som];

% performing linear initialization of nodes
display(['initialization ' num2str(ny_som) 'x' num2str(nx_som) ' SOM'])
sMap=som_lininit(data,'msize',msize,'hexa','sheet');

% training SOM
display('training')
[sM,sT] = som_batchtrain(sMap,data,'ep','hexa','sheet','radius',[3 1],'trainlen',200); 

% calculating quantization and topological error
[QE(k),TE(k)]=som_quality(sM,data)
% calulating hits (frequencies) of occurences of each pattern
hi=som_hits(sM,data);
hi=100*hi/sum(hi);

% plotting the sMap and SOM
index=NaN; 
index=[1:en];
index=reshape(index,ny_som,nx_som);
index=index';
index=reshape(index,1,nx_som*ny_som);

% plot the final SOM patterns
figure;
for i=1:en
subplot(ny_som,nx_som,i);
Y1=reshape(sM.codebook(index(i),:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);
title([num2str(index(i)) ', f=' num2str(hi(index(i)),'%2.1f')])
axis off
end

% alternative way of plotting SOM with faces 
numPerLine = nx_som;
ShowLine = ny_som;

Y = zeros(faceH*ShowLine,faceW*numPerLine);
for i=0:ShowLine-1
   for j=0:numPerLine-1
     Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW) = reshape(sM.codebook(index(i*numPerLine+j+1),:),[faceH,faceW]);
   end
end

figure;
imagesc(Y);
colormap(gray);
axis off
title([num2str(ny_som) 'x' num2str(nx_som) ' SOM'])
end

% plotting QE and TE
figure; 
plotyy([1:k],QE,[1:k],TE);
[hAx,hLine1,hLine2] = plotyy([1:k],QE,[1:k],TE);
xlim(hAx(1),[1 k]);
xlim(hAx(2),[1 k]);
set(hAx(1),'XTick',[1:1:6],'XTicklabel',{'3x3','4x3','5x3','4x4','5x4','5x5'});
set(hAx(2),'XTick',[1:1:6],'XTicklabel',{'3x3','4x3','5x3','4x4','5x4','5x5'});
ylabel(hAx(1),'QE') % left y-axis
ylabel(hAx(2),'TE') % right y-axis

% according to this plot one can decide to go for 4x4 SOM as an optimal size
% but in general it depends on the user ho much level of detail to keep:
% larger SOM sizes retain more details

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDITIONAL (optional) material 
% the following algorithm determines optimal SOM size, however it's usualy very large SOM
sD = som_data_struct(data); 
sM = som_make(sD); 
ny_som = sM.topol.msize(1);
nx_som = sM.topol.msize(2);
en=ny_som*nx_som;

Bmus = som_bmus(sM,sD);

% calculating hits (frequencies) of occurences of each pattern
hi=som_hits(sM,data);
hi=100*hi/sum(hi);

% plotting the sMap and SOM
index=[1:en];
index=reshape(index,ny_som,nx_som);
index=index';
index=reshape(index,1,nx_som*ny_som);

% plotting SOM with faces 
numPerLine = nx_som;
ShowLine = ny_som;

Y = zeros(faceH*ShowLine,faceW*numPerLine);
for i=0:ShowLine-1
   for j=0:numPerLine-1
     Y(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW) = reshape(sM.codebook(index(i*numPerLine+j+1),:),[faceH,faceW]);
   end
end

figure;
imagesc(Y);
colormap(gray);
axis off
title([num2str(ny_som) 'x' num2str(nx_som) ' SOM'])


