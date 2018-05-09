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

% 1) In the space of first few eigenvectors (decide yourself how many) perform the hierarchical clustering and find the optimal number of clusters, i.e. decide the optimal number from the dendrogram. 

variance=eigenvalues./sum(eigenvalues);

figure; 
plot([1:length(variance)],variance,'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');
% from this figure I decided to keep first two modes for further analysis

% clustering in the space of first 3 modes
PCtot=PCs(:,1:2);

% find the dendogram
Zn = linkage(PCtot,'ward','euclidean');

figure;
dendrogram(Zn);
% from dendrogram I decided to go for 5 clusters (also could chose 3 or 4)

Nc=5;
cn = cluster(Zn,'maxclust',Nc);

figure;
scatter(PCs(:,1),PCs(:,2),30,cn,'filled')
xlabel('PC1')
ylabel('PC2');

%2) Find the centre points (arithmetic mean) of each cluster in the space of eigenvectors and reconstruct the face for each of this centre points. Compare this reconstructed face with the mean reconstructed face from each cluster. [Hint: for the latter, first reconstruct all the images based on the first few modes and then calculate the mean image for each cluster -> see the MATLAB code below]

% calculate the centre points of each cluster in the PC space 
for i=1:Nc
[ind dummy]=find(cn == i);
PCs_center(i,:)=mean(PCs(ind,1:2));
end

% reconstructing the face for those centre points 
for j=1:Nc
data_rec1_mean(j,:)=PCs_center(j,1).*eigenvectors(:,1)+PCs_center(j,2).*eigenvectors(:,2); 
end 

% plot the reconstructed faces of those centre points
figure; 
for i=1:Nc
subplot(2,3,i)
Y1=reshape(data_rec1_mean(i,:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);
title(['Rec. face cluster center ',num2str(i)]);
end

% reconstructing all the faces from the first two PCs
for j=1:size(PCs,1);
data_rec2(j,:)=PCs(j,1).*eigenvectors(:,1)+PCs(j,2).*eigenvectors(:,2); 
end 

% calculate the mean reconstructed face for each cluster 
for i=1:Nc
[ind dummy]=find(cn == i);
data_rec2_mean(i,:)=mean(data_rec2(ind,:));
end

% plot the mean reconstructed face for each cluster
figure; 
for i=1:Nc
subplot(2,3,i)
Y1=reshape(data_rec2_mean(i,:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);
title(['Mean rec. face ',num2str(i)]);
end

% plot the difference b/w data_rec1_mean and data_rec2_mean
delta=data_rec1_mean-data_rec2_mean;
figure; 
for i=1:Nc
subplot(2,3,i)
Y1=reshape(delta(i,:),[faceH,faceW]);
imagesc(Y1);
colormap(gray);
colorbar
title(['Difference for cluster ',num2str(i)]);
end

% can see that the differences above are negligible (difference are of the order 10^(-16)), 
% so either method gives the same results in this example


%3) Using the same optimal number of clusters perform the K-means clustering on the same dataset (in PC space). How much do the K-means centre points differ from the centre points derived from the hierarchical clustering? 

% K-means clustering
K=Nc;
opts = statset('Display','final');
[idx,K_center] = kmeans(PCs(:,1:2),K,'Replicates',20,'Options',opts);  % here choosing 20 runs

% plotting the centre points from hierarchical clustering and from the K-means
figure;
subplot(1,2,1)
scatter(PCs(:,1),PCs(:,2),30,cn); hold on 
plot(PCs_center(:,1),PCs_center(:,2),'b*'); 
xlabel('PC1')
ylabel('PC2');
legend('clustered data','cluster centers');

subplot(1,2,2)
scatter(PCs(:,1),PCs(:,2),30,idx); hold on
plot(K_center(:,1),K_center(:,2),'r*');
xlabel('PC1')
ylabel('PC2');
legend('K-means clustered data','K-mean centers');

% we can see that the clusters a quite different (differently colored points);
% lets calculate the minimum distances between each cluster center and K-mean centers

% Euclidean distance between two centers:
for j=1:Nc
  for k=1:K
D(k)=sqrt(sum((PCs_center(j,:)-K_center(k,:)).^2));
  end
[dist(j) index(j)]=min(D);
clear D
end

% for each cluster index plot its matching (minimum distance)
% K-means cluster index, and then
% plot the distances between the two cluster centers 
figure;
subplot(2,1,1);
plot([1:Nc],index,'bo-');
xlabel('cluster index');
ylabel('K-means cluster index');

subplot(2,1,2);
plot([1:Nc],dist,'ro-');
xlabel('cluster index');
ylabel('distance to K-means point');











