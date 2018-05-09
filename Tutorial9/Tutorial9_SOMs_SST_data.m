%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 9 (31 Oct 2017)
% SOMs on SST data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2:
% This is the same data as from Tutorial 4 on PCA:
% Gridded monthly sea surface temperature data for Tropical Pacific
% from ERA Interim reanalysis
% Period: Jan 1979 to Jun 2015

% NOTE: in order to run this script you need to have
% SOM Toolbox (downloaded from the website:
% http://www.cis.hut.fi/projects/somtoolbox/)
% included in your Matlab path 
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
xlim([1979 2015.6]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOM algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initilizing SOM -> this takes several minutes depending on the speed of your comp
% chose size of the map for SOM 
% I tried only 3x4 SOM
ny_som=3; nx_som=4;
en=ny_som*nx_som;

data=double(data);

msize=[ny_som nx_som];
% performing linear initialization of nodes
display('initialization')
sMap=som_lininit(data,'msize',msize,'hexa','sheet');

% training SOM
display('training')
[sM,sT] = som_batchtrain(sMap,data,'ep','hexa','sheet','radius',[3 1],'trainlen',200); 

% calulating quantization error
[q,t]=som_quality(sM,data)

% calulating hits (frequencies) of occurences of each pattern, for each seasn
hi=som_hits(sM,data);
hi=100*hi/sum(hi);

% plotting the sMap and SOM
index=[1:en];
index=reshape(index,ny_som,nx_som);
index=index';
index=reshape(index,1,nx_som*ny_som);

figure;
for i=1:en
subplot(ny_som,nx_som,i);
var01=sM.codebook(index(i),:);
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
if i == en
colorbar;
end
caxis([min(min(sM.codebook)) max(max(sM.codebook))])
title([num2str(index(i)) ', f=' num2str(hi(index(i)),'%2.1f')])
end

% find the best-matching units from the map for the given vectors
bmus=som_bmus(sM,data);

% plot SOM node index in time
figure;
plot(time,bmus,'bo-');
xlabel('Time (yr)');
ylabel('SOM node index');

% extend bmus with NaN so that is starts with Jan 1979 to Dec 2017
bmusn(1:39*12)=NaN;
bmusn(2:end-7)=bmus;

% plot bmus in 2D
x1=[1:12];
y1=[1979:2017];
figure;
bmus_plot=reshape(bmusn,12,39);
h=imagesc(y1,x1,bmus_plot);
colormap jet
colorbar
set(h,'alphadata',~isnan(bmus_plot));
xlabel('Time (yr)');
ylabel('Month');


