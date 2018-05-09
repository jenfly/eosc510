%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 8 (24 Oct 2017) 
% Clustering of synthetic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1
% apply hierarchical clustering (Ward's method) to
% synthetic data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the data;
load fisheriris;

% create the data for this example 
meas(101:150,3)=meas(1:50,4);
meas=meas(1:150,1:3);

% same data but randomized:
indi=randperm(150);
meas_rand=meas(indi,:);

% find links between each data pair and plot the dendrogram:
Z = linkage(meas,'ward','euclidean');

figure;
dendrogram(Z)
title('dendrogram');

% cluster the data by choosing the number of clusters k
k=3;
c = cluster(Z,'maxclust',k);

% plot original data and the clusters
figure;
subplot(2,1,1);
imagesc(meas_rand');
colormap jet
colorbar
title('original data');
xlabel('points');
ylabel('measured variables');

subplot(2,1,2);
imagesc(c')
colormap jet
colorbar
title('Clusters');
xlabel('points');

% plot original data in 3D
figure;
scatter3(meas(:,1),meas(:,2),meas(:,3),30,'filled');
xlabel('x');
ylabel('y');
zlabel('z');

% plot clustered data in 3D
% each cluster in different color
figure;
scatter3(meas(:,1),meas(:,2),meas(:,3),30,c,'filled')
xlabel('x');
ylabel('y');
zlabel('z');
title('Clusters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same example but using K-means clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K-means clustering
% select K (number of clusters):
K=3;
opts = statset('Display','final');
[idx,ctrs] = kmeans(meas,K,'Replicates',10,'Options',opts);

% plot original data and the clusters
figure;
subplot(2,1,1);
imagesc(meas');
colormap jet
colorbar
title('original data');
xlabel('points');
ylabel('measured variables');

subplot(2,1,2);
imagesc(idx')
colormap jet
colorbar
title('K-means clusters');
xlabel('points');

% plot clusters in 3D and the cluster centers (ctrs)
figure;
scatter3(meas(:,1),meas(:,2),meas(:,3),30,idx,'filled'); 
hold on
plot3(ctrs(:,1),ctrs(:,2),ctrs(:,3),'r*'); 
xlabel('x');
ylabel('y');
zlabel('z');
title('K-means clusters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2
% apply hierarchical clustering (Ward's method) to
% a synthetic data in combination with PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating the input data
k=2*pi/100;
omega=2*pi/50;
x=[1:100];
t=[1:200];

for i=1:200
y(i,:)=sin(k*x-omega*t(i));
end

% sine pattern w/ amplitude of 1
y1=y(1,:);

% cosine pattern w/ amplitude of 0.5
y2=0.5*y(37,:);

% step pattern w/ amplitude of 0.8
y3(1:50)=-0.8;
y3(51:100)=0.8;

% sawtooth pattern w/ amplitude of 1
y4(1:25)=-2*x(1:25)./25+1;
y4(26:50)=2*x(26:50)./25-3;
y4(51:75)=-2*x(51:75)./25+5;
y4(76:100)=2*x(76:100)./25-7;

% plot the patterns in the data 
% before the noise is added
figure;
subplot(2,2,1)
plot(y1);
ylim([-1 1]);
title('(a)');
xlabel('x');
ylabel('y');

subplot(2,2,2)
plot(y3);
ylim([-1 1]);
title('(b)');
xlabel('x');
ylabel('y');

subplot(2,2,3)
plot(y4);
ylim([-1 1]);
title('(c)');
xlabel('x');
ylabel('y');

subplot(2,2,4)
plot(y2);
ylim([-1 1]);
title('(d)');
xlabel('x');
ylabel('y');

% create timeseries where from t=1:50 signal y1 occurs,
% then for the next 50 steps signal y3 occurs,
% then for the next 50 steps signal y4 occurs,
% and then for the last 50 steps signal y2 occurs
ynew(1:50,:)=repmat(y1,50,1);
ynew(51:100,:)=repmat(y3,50,1);
ynew(101:150,:)=repmat(y4,50,1);
ynew(151:200,:)=repmat(y2,50,1);

% add noise to the data
load 'noise.mat';
ynew1=ynew+2*noise;

ynew=[ynew1; ynew1];

% plot few snapshots of the data
labels={'(a)','(b)','(c)','(d)'};
figure;
for i=1:4
subplot(2,2,i)
plot(ynew(i*50-1,:));
ylim([-2 2]);
title(labels(i));
xlabel('x');
ylabel('y');
end

% PCA:
[eigenvectors,PCs,eigenvalues]=princomp(ynew);

% contribution of each mode to total variance:
variance=eigenvalues./sum(eigenvalues);

figure; 
plot([1:100],variance,'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

% plot first three modes (eigenvectors and PCs)
figure; 
for i=1:3
subplot(3,2,2*i-1)
plot(eigenvectors(:,i));
%ylim([-1 1]);
xlabel('x');
title(['eigenvector',num2str(i)]);

subplot(3,2,2*i)
plot(PCs(:,i));
xlabel('time');
title(['PC',num2str(i)]);
end

% plot the data in the space of first three eigenvectors
figure;
scatter3(PCs(:,1),PCs(:,2),PCs(:,3));
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('data in the space of first 3 eigenvectors');

% put first 3 PCs into one matrix: 
PCtot=PCs(:,1:3);

% perfrom clustering on this new matrix PCtot
% find the dendogram
Zn = linkage(PCtot,'ward','euclidean');

figure;
dendrogram(Zn);

% select number of clusters
k=4;
cn = cluster(Zn,'maxclust',k);

% plot the clustered data in the space of first three eigenvectors
figure;
scatter3(PCs(:,1),PCs(:,2),PCs(:,3),16,cn)
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('clustered data');

% calculate mean signal (spatial pattern) for each cluster
for i=1:4
[ind dummy]=find(cn == i);
y_mean_model(i,:)=mean(ynew(ind,:));
end

% calulate mean signal from original data (where we know that there are 4 clusters)
y_mean_obs(1,:)=mean(ynew([1:50 201:250],:));
y_mean_obs(2,:)=mean(ynew([51:100 251:300],:));
y_mean_obs(3,:)=mean(ynew([101:150 301:350],:));
y_mean_obs(4,:)=mean(ynew([151:200 351:400],:));

% get the order right, so that the order of modeled mean signals corresponds 
% to the order of original mean signals: 
for i=1:4
r=corrcoef([y_mean_obs(i,:)' y_mean_model']);
[maxi index(i)]=max(r(1,2:end));
end

% compare original mean signal with modeled one 
figure;
for i=1:4
subplot(2,2,i)
plot([1:100],y_mean_obs(i,:),'b-',[1:100],y_mean_model(index(i),:),'r-');
legend('original mean','reconstructed mean');
ylim([-2 2])
end


