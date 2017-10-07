%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 3 (21 Sep 2016) with solutions 
% PCA on human face data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This dataset contains grayscale images of 38 individuals and around 64 near frontal images under different illuminations per individual. The faces are processed so that they are cropped and centered, with 32 x 32 pixels each.

load YaleB_32x32.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scale the features (pixel values) to [0,1]
maxValue = max(max(fea));
fea = fea/maxValue;
%===========================================

% plot 100 first images all together on one figure
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
title('original images');

% example of a plot (just the first face) 
figure;
Y1=reshape(fea(1,:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);
title('the first image in the set');

%1) Find the eigenfaces
% apply PCA:
%option 1:
[eigenvectors,PCs,eigenvalues]=princomp(fea);


% option 2: apply PCA mannualy
x=fea';
%subtract mean
x=bsxfun(@minus, x, mean(x,2));
% calculate covariance
s = cov(x');
% obtain eigenvalue & eigenvector
[V,D] = eig(s);
eigval = diag(D);
% sort eigenvalues in descending order
eigval = eigval(end:-1:1);
V = fliplr(V);
% get principal components
a=V'*x;

%2) How many modes are enough too represent the images of human faces?
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
% from this figure I decided to keep 7 first modes 
% any answer from 2 to 8 modes is fine

% plot the first seven modes (eigenvectors and PCs)
figure; 
for i=1:7
subplot(2,7,i)
Y1=reshape(eigenvectors(:,i),[faceH,faceW]);

%ylim([-1 1]);
imagesc(Y1);
colormap(gray);
title(['e_',num2str(i)]);

subplot(2,7,i+7)
plot(PCs(:,i));
xlabel('index');
title(['PC',num2str(i)]);
xlim([1 length(gnd)])
end

% reconstructed faces from the first k modes
k=7;
invV=inv(V');
y_rec=(invV(:,1:k)*a(1:k,:))';

% plot 100 first reconstructed faces all togeter on one figure
Yrec = zeros(faceH*ShowLine,faceW*numPerLine);
for i=0:ShowLine-1
   for j=0:numPerLine-1
     Yrec(i*faceH+1:(i+1)*faceH,j*faceW+1:(j+1)*faceW) = reshape(y_rec(i*numPerLine+j+1,:),[faceH,faceW]);
   end
end

figure;
imagesc(Yrec);
colormap(gray);
title('reconstructed images');

%3) What is the most 'generic' face in the set of images?
% finding the 'generic' face, i.e. the one with minimum RMSE
err=fea-y_rec;
RMSE=(sum((err.^2)')./(32*32)).^0.5;
[minRMSE indmin]=min(RMSE);

% plot the most generic face
figure;
Y1=reshape(fea(indmin,:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);
title('the most generic face')

% we can see that the most generic face is almost totaly shaded, 
% so lets plot, for example, 14 faces with the smallest RMSE (sorted in ascending order of RMSE)
[sortRMSE indsort]=sort(RMSE);
figure; 
for i=1:14
subplot(2,7,i)
Y1=reshape(fea(indsort(i),:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);
title(['Face ',num2str(i)]);
end


