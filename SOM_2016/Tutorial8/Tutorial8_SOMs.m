%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 8 (26 Oct 2016)
% Application of SOM on a synthetic data 
% Example 1:
% compare SOM performance with the clustering algorithm
% and PCA
% NOTE: in order to run this script you need to have
% SOM Toolbox installed from the website:
% http://www.cis.hut.fi/projects/somtoolbox/
% in your Maltab path include the downloaded folder with som toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../somtoolbox')

load 'data.mat';

% plot the data (all threee variables) in time
figure;
imagesc(data');
xlabel('time');
colorbar
title('3-D data in time');

% basic features in the data (mean features)
data_basic(1,:)=mean(data(1:50,:));
data_basic(2,:)=mean(data(51:100,:));
data_basic(3,:)=mean(data(101:150,:));
data_basic(4,:)=mean(data(151:200,:));

% plot those mean features
figure; 
for i=1:4
subplot(2,2,i)
plot([1:3],data_basic(i,:),'bo-')
ylim([-8 10]);
title(['mean pattern' num2str(i)])
end

% plot the data in 3-D 
figure;
scatter3(data(:,1),data(:,2),data(:,3));
xlabel('x1');
ylabel('x2');
zlabel('x3');
title('data in 3-D')

% plot bunch of data
figure;
for i=1:10
subplot(4,10,i);
plot([1:3],data(i,:),'ko-');
ylim([-8 10]);
subplot(4,10,i+10);
plot([1:3],data(i+50,:),'ko-');
ylim([-8 10]);
subplot(4,10,i+20);
plot([1:3],data(i+100,:),'ko-');
ylim([-8 10]);
subplot(4,10,i+30);
plot([1:3],data(i+150,:),'ko-');
ylim([-8 10]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply clustering algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find the dendogram
Zn = linkage(data,'ward','euclidean');

figure;
dendrogram(Zn);
title('dendogram');

% find four clusters
cn = cluster(Zn,'maxclust',4);

C=[1 0 1; 0 1 0; 0 0 1; 0 1 1];
Cn=C(cn,:);

% plot the data with colored clusters
figure;
scatter3(data(:,1),data(:,2),data(:,3),16,Cn)
title('3-D data with clusters');

% calculate mean signal for each cluster
for i=1:4
[ind dummy]=find(cn == i);
data_cluster(i,:)=mean(data(ind,:));
end

% get the order right, so that the model signal order corresponds to the original signal order
for i=1:4
r=corrcoef([data_basic(i,:)' data_cluster']);
[maxi index(i)]=max(r(1,2:end));
end

% compare original mean signal with modeled (clustered) one 
figure;
for i=1:4
subplot(2,2,i)
plot([1:3],data_basic(i,:),'bo-',[1:3],data_cluster(index(i),:),'ro-');
legend('mean original','clustered');
ylim([-8 10])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA:
[eigenvectors,PCs,eigenvalues]=princomp(data);

% contribution of each mode to total variance:
variance=eigenvalues./sum(eigenvalues);

figure; 
plot([1:3],variance,'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

% plot first three modes (eigenvectors and PCs)
figure; 
for i=1:3
subplot(3,2,2*i-1)
plot([1:3],eigenvectors(:,i),'bo-');
title(['eigenvector',num2str(i)]);

subplot(3,2,2*i)
plot(PCs(:,i));
xlabel('time');
title(['PC',num2str(i)]);
end

% plot the data with eignevectors
% the eignevector coincide with the coordinates in the scatter plot
f=[0:10];
g(1:11)=0;
figure;
scatter3(data(:,1),data(:,2),data(:,3),16,Cn)
hold on
line(f,g,g); hold on
line(g,f,g); hold on
line(g,g,f); 
xlabel('e1');
ylabel('e2');
zlabel('e3');
title('data with eigenvectors);

%%%%%%%%%%%%%%%%%%%%%%%%%
% apply SOM -> for this you need SOM Toolbox
%%%%%%%%%%%%%%%%%%%%%%%

% initilizing SOM
% chose size of the map for SOM 
%ny_som=2; nx_som=2;
%ny_som=2; nx_som=3;
ny_som=3; nx_som=2;
en=ny_som*nx_som;

data=double(data);

msize=[ny_som nx_som];
% performing linear initialization of nodes (i.e. the nodes are distibuted in the space of 
% firt two eignevectors from PCA
display('initialization')
sMap=som_lininit(data,'msize',msize,'hexa','sheet');

% training SOM
display('training')
[sM,sT] = som_batchtrain(sMap,data,'ep','hexa','sheet','radius',[2 1],'trainlen',200); 

% calulating quantization error
[q,t]=som_quality(sM,data)

% calulating hits (frequencies) of occurences of each pattern, for each seasn
hi=som_hits(sM,data);
hi=100*hi/sum(hi);

% plotting the sMap (initialized map) and final SOM
% here get the order right for the plot
index=[1:en];
index=reshape(index,ny_som,nx_som);
index=index';
index=reshape(index,1,nx_som*ny_som);

% plot inital pattterns 
figure;
for i=1:en
subplot(ny_som,nx_som,i);
plot([1:3],sMap.codebook(index(i),:),'bo-');
xlabel('time');
title(['initial node ' num2str(index(i))])
ylim([-8 10])
end

% plot final SOM patterns
figure;
for i=1:en
subplot(ny_som,nx_som,i);
plot([1:3],sM.codebook(index(i),:),'bo-');
xlabel('time');
title(['SOM node ' num2str(index(i)) ', f=' num2str(hi(index(i)),'%2.1f')])
ylim([-8 10])
end

% find the order of the SOM nodes in time 
% and them plot SOM patterns in time
bmus=som_bmus(sM,data);
figure; 
imagesc(sM.codebook(bmus,:)');
colorbar;
xlabel('time');
title('SOM patterns in time');

% plot inital grid with the data in 3-D
figure;
scatter3(data(:,1),data(:,2),data(:,3),16,Cn)
hold on
line(f,g,g); hold on
line(g,f,g); hold on
line(g,g,f); 
xlabel('e1');
ylabel('e2');
zlabel('e3');
scatter3(sMap.codebook(:,1),sMap.codebook(:,2),sMap.codebook(:,3),30,[1 0 0],'filled')
title('Initial grid');

% plot final SOM points with the data in 3-D
figure;
scatter3(data(:,1),data(:,2),data(:,3),16,Cn)
hold on
line(f,g,g); hold on
line(g,f,g); hold on
line(g,g,f); 
xlabel('e1');
ylabel('e2');
zlabel('e3');
scatter3(sM.codebook(:,1),sM.codebook(:,2),sM.codebook(:,3),30,[1 0 0],'filled')
title('Final SOM');





