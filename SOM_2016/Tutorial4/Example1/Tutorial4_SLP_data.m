%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 4 (28 Sep 2015)
% PCA on real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1:
% Gridded daily sea-level pressure data for North Pacific domain
% from CFSR climate reanalysis
% Period: 1 Jan 1979 to 31 Dec 2010
% the .mat file is organized so that rows are grid cells (representing SLP
% value (Pa) for given latitude and logitude)
% while columns are days (1979-2010)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first loading the data and plotting it
% bunch of files and commands here are just for the plots 

load 'coordinates_for_plots.mat';
load 'incrop.mat';

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
filename='SLP_CFSR_daily_200km_cropped_1979_2010.mat';  
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

figure; 
subplot(2,1,1)
plot([1:length(incrop)],variance,'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

subplot(2,1,2)
plot([1:20],variance(1:20),'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

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

% plot data in space of the two first two PCs
figure;
plot(PCs(indicesw,1),PCs(indicesw,2),'b.',PCs(indicess,1),PCs(indicess,2),'r.');
xlabel('PC1');
ylabel('PC2');
legend('winter','summer')
grid on
xlim([-800 800])
ylim([-800 800])





