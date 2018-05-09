%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 8 (24 Oct 2017)
% clustering on real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 3:
% This is the same data as from Tutorial 4 on PCA:
% Gridded monthly sea surface temperature data for Tropical Pacific
% from ERA Interim reanalysis
% Period: Jan 1979 to Jun 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load SST anomalies (the processed data from Tutorial 4)
load SST_anomalies_data.mat;

% these variables are needed for plots:
tt=[135 29 462]; % number of longitudes, latitudes, and points in time
x=[159.75:0.75:260.25]';  % longitudes
y=[10.50:-0.75:-10.50]';  % latitudes
time=[1979+2/12:1/12:2017+5/12];

% PCA:
[eigenvectors,PCs,eigenvalues]=princomp(data);

% contribution of each mode to total variance:
variance=eigenvalues./sum(eigenvalues);

% plot first 4 modes (eigenvectors and PCs)
figure; 
for i=1:4
subplot(4,2,2*i-1)
var01=eigenvectors(:,i);
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
colorbar;
title(['eigenvector var=',num2str(variance(i),'%2.2f')]);

subplot(4,2,2*i)
plot(time,PCs(:,i),'b-')
xlabel('time');
title(['PC',num2str(i)]);
xlim([1979 2017.6]);
end

% plot the data in the space of first three eigenvectors
figure;
scatter3(PCs(:,1),PCs(:,2),PCs(:,3),30,'filled');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('original data');

% create a matrix with only first 3 PCs
PCtot=PCs(:,1:3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLUSTERING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the dendogram on this PCtot matrix
Zn = linkage(PCtot,'ward','euclidean');

figure;
dendrogram(Zn);

% find k clusters
k=6;
cn = cluster(Zn,'maxclust',k);

% plot the clusters in the space of first 3 eigenvectors
figure;
scatter3(PCs(:,1),PCs(:,2),PCs(:,3),30,cn,'filled')
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('clustered data');

% plot the clusters in the space of first 2 eigenvectors
figure;
scatter(PCs(:,1),PCs(:,2),30,cn,'filled')
xlabel('PC1');
ylabel('PC2');
title('clustered data');

% reconstructing all data from the first 3 modes 
for j=1:length(data(:,1))
data_rec(j,:)=PCs(j,1).*eigenvectors(:,1)+PCs(j,2).*eigenvectors(:,2)+PCs(j,3).*eigenvectors(:,3);
end

% calculate the mean (spatial) pattern of each cluster
for i=1:k
[ind dummy]=find(cn == i);
data_mean_model(i,:)=mean(data_rec(ind,:));
end

% plot the mean pattern of each cluster
figure; 
for i=1:k
subplot(2,3,i)
var01=data_mean_model(i,:);
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
if i==k
colorbar;
end
caxis([-1.5 3])
title(['cluster ',num2str(i)]);
end

% plot the time-series of clusters
figure;
plot(time,cn,'k-'); hold on
scatter(time,cn,30,cn,'filled'); hold on
xlabel('time');
ylabel('Cluster');
xlim([1979 2017.6]);
ylim([0 k+1]);

