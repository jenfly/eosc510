%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 9 (2 Nov 2016)
% SOMs on real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1:
% This is the same data as from Tutorial 4 on PCA:
% Gridded daily sea-level pressure data for North Pacific domain
% from CFSR climate reanalysis
% Period: 1 Jan 1979 to 31 Dec 2010
% The first part of the script is a copy from Tuturial4_SLP_data.m
% NOTE: in order to run this script you need to have
% SOM Toolbox installed from the website:
% http://www.cis.hut.fi/projects/somtoolbox/
% in your Maltab path include the downloaded folder with som toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first loading the data and plotting it
% bunch of files and commands here are just for the plots 
% create your path here for the directory where the data from Tutorial 4 is (these paths are for my computer):
load '/shared/EOSC510/Tutorials/Tutorial4/Example1/coordinates_for_plots.mat';
load '/shared/EOSC510/Tutorials/Tutorial4/Example1/incrop.mat';

helpers=[1403 1457 1459 1513];
ind_bad=1458;

[IN1 ON1]=inpolygon(x_coast_NARR,y_coast_NARR,x_domp,y_domp);
xx_coast_NARR=x_coast_NARR;
yy_coast_NARR=y_coast_NARR;
xx_coast_NARR(IN1 == 0)=NaN;
yy_coast_NARR(IN1 == 0)=NaN;

X1=reshape(x_NARR,Nx,Ny);
Y1=reshape(y_NARR,Nx,Ny);

X1b=X1(1:5:end,1:5:end);
Y1b=Y1(1:5:end,1:5:end);

ylow=-1e6;
yhigh=1e6;
xlow=-10e5;
xhigh=18e5;

%%%%%%%%%%%%%%%%%%%%
% loading CFSR SLP file
%%%%%%%%%%%%%%%%%%
filename='/shared/EOSC510/Tutorials/Tutorial4/Example1/SLP_CFSR_daily_200km_cropped_1979_2010.mat';  
load(filename);
SLP(ind_bad,:)=mean(SLP(helpers,:));
var0=SLP;

% finding indices for winter and summer 
indices=NaN;
ni = datenum(1979,1,1);
index=[1:11688]+ni-1;

% winter indices
for year=1979:2009
d1=datenum(year,12,1);
d2=datenum(year+1,3,1);

[dummy ind1]= find(index == d1);
[dummy ind2]= find(index == d2);

index1=[ind1:ind2];
indices=[indices index1];
end
indicesw=indices(2:end);

% summer indices
indices=NaN;
for year=1979:2010
d1=datenum(year,6,1);
d2=datenum(year,9,1);

[dummy ind1]= find(index == d1);
[dummy ind2]= find(index == d2);

index1=[ind1:ind2];
indices=[indices index1];
end
indicess=indices(2:end);

VAR0=var0(incrop,:)./100;  % setting SLP to hPa units
var01(1:Nx*Ny)=NaN;

% plot SLP over the whole domain
% for first 30 days in Jan 1979

figure;
for ii=1:30
var01(incrop)=VAR0(:,ii);
var_plot=reshape(var01,Nx,Ny);
subplot(5,6,ii)
h=imagesc(x_NARR,y_NARR,var_plot');
if ii==30
colorbar('SouthOutside' );
end
colormap jet
caxis([min(min(VAR0(:,1:30))) max(max(VAR0(:,1:30)))])
axis xy  
set(h,'alphadata',~isnan(var_plot'))  
hold on
plot(xx_coast_NARR,yy_coast_NARR,'k-','LineWidth',1);
hold on
plot(x_dom,y_dom,'-.','Color',[0.6 0.6 0.6],'LineWidth',1);
xlim([min(x_dom),max(x_dom)]);
ylim([min(y_dom),max(y_dom)]); 
title(['Jan ',num2str(ii),' 1979'])
axis off 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparing the data for PCA analysis
data=VAR0';

% PCA:
[eigenvectors,PCs,eigenvalues]=princomp(data);

% contribution of each mode to total variance:
variance=eigenvalues./sum(eigenvalues);

% set the time series (in years)
time=[1979:1/365.25:2011-1/365.25];

% separate indices for winter and summer (just so that in the plots
% they are better visualized
nn=size(PCs);
PCs_winter(1:nn(1),1:nn(2))=NaN;
PCs_summer(1:nn(1),1:nn(2))=NaN;
PCs_winter(indicesw,:)=PCs(indicesw,:);
PCs_summer(indicess,:)=PCs(indicess,:);

% plot first 4 modes (eigenvectors and PCs)
figure; 
for i=1:4
subplot(4,2,2*i-1)
var01(incrop)=eigenvectors(:,i);
var_plot=reshape(var01,Nx,Ny);
h=imagesc(x_NARR,y_NARR,var_plot');
colormap jet
colorbar;
%caxis([min(min(eigenvectors)) max(max(eigenvectors))])
axis xy  
set(h,'alphadata',~isnan(var_plot'))  
hold on
plot(xx_coast_NARR,yy_coast_NARR,'k-','LineWidth',1);
hold on
plot(x_dom,y_dom,'-.','Color',[0.6 0.6 0.6],'LineWidth',1);
xlim([min(x_dom),max(x_dom)]);
ylim([min(y_dom),max(y_dom)]); 
title(['eigenvector var=',num2str(variance(i),'%2.2f')]);
axis off 

subplot(4,2,2*i)
plot(time,PCs(:,i),'y-',time,PCs_winter(:,i),'b-',time,PCs_summer(:,i),'r-')
xlabel('time');
title(['PC',num2str(i)]);
xlim([1979 2011]);
legend('all','winter','summer');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOM algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var01(1:Nx*Ny)=NaN;
% initilizing SOM -> this takes several minutes depending on the speed of your comp
% chose size of the map for SOM  
% the code is made for 4x5, 5x7, 3x4 
ny_som=4; nx_som=5;
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
var01(incrop)=sM.codebook(index(i),:);
var_plot=reshape(var01,Nx,Ny);
h=imagesc(x_NARR,y_NARR,var_plot');
colormap jet
if i == en
colorbar;
end
caxis([min(min(sM.codebook)) max(max(sM.codebook))])
axis xy  
set(h,'alphadata',~isnan(var_plot'))  
hold on
plot(xx_coast_NARR,yy_coast_NARR,'k-','LineWidth',1);
hold on
plot(x_dom,y_dom,'-.','Color',[0.6 0.6 0.6],'LineWidth',1);
xlim([min(x_dom),max(x_dom)]);
ylim([min(y_dom),max(y_dom)]); 
title([num2str(index(i)) ', f=' num2str(hi(index(i)),'%2.1f')])
axis off 
end

% find the best-matching units from the map for the given vectors
bmus=som_bmus(sM,data);

%getting rid of the leap years
indexN=[1:11688]+ni-1;
br=0;
for year=1980:4:2010
d1=datenum(year,2,29);
br=br+1;
[dummy ind1(br)]= find(indexN == d1);
end

bmus_new=bmus;
bmus_new(ind1)=NaN;

% this is bmus w/o 29 Feb for leap years
bmusn=bmus_new(~isnan(bmus_new));

% plot bmus in 2D
x1=[1:365];
y1=[1979:2010];
figure;
imagesc(y1,x1,reshape(bmusn,365,32));
colormap jet
colorbar;
xlabel('Time (yr)');
ylabel('Day of year');

% plot colored nodes in SOM
if en == 20
imi=[1:4; 5:8; 9:12; 13:16; 17:20]; 
elseif en == 35
imi=[1:5; 6:10; 11:15; 16:20; 21:25; 26:30; 31:35];
else
imi=[1:3; 4:6; 7:9; 10:12];
end

figure;
imagesc(imi');
colormap jet
colorbar;
axis off

