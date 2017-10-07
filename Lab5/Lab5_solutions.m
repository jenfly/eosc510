%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lab 5 (5 Oct 2017)
% CCA on real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code loads SST anomalies data from last tutorial (saved in SST_anomalies_data.mat)
% then load .nc file of gridded monthly temperature data for BC
% Period: Jan 1979 to July 2017
% all the data pre-processing is done below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load SST anomalies (the processed data from last Tutorial)
load SST_anomalies_data.mat;
data_SST=data;
% these are needed for plots:
tt_SST=[135 29 462];  % contains number of longitude, latitude and time points for SST data
x_SST=[159.75:0.75:260.25]';  % longitudes for SST data
y_SST=[10.50:-0.75:-10.50]';  % latitudes for SST data
% put this data into variable xdata
xdata=data_SST;

months={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};

% loading BC temp data
% this is nc formated file (binary, structured file)
filename='T_monthly_Canada_1979_2017_ERA_Interim.nc';

% the following command will display the file content
ncdisp(filename);

% short way of reading the .nc file:
var=ncread(filename,'t2m');
x=ncread(filename,'longitude');
y=ncread(filename,'latitude');
time=ncread(filename,'time');

tt=size(var);
var=double(var);
var=var-273.15; % setting units to degrees Celsius  

% var is 3-D matrix, here I am reshaping it into 2-D
% where rows are grid points, columns are time points
for kj=1:tt(3)
var1=squeeze(var(:,:,kj));
vard(1:tt(1)*tt(2),kj)=reshape(var1,tt(1)*tt(2),1);
end

% converting integer values to double precision 
vard=double(vard);

% lets plot temp over the whole domain (BC) for Jan 1979
var_plot=reshape(vard(:,1),tt(1),tt(2));
load coastline.mat
figure;
h=imagesc(x,y,var_plot');
colormap jet
colorbar
axis xy
hold on
plot(long,lat,'k');
plot(long+360.,lat,'k');
ylim([45 60]); 
xlim([220 245])

% plot temp for all months in 1988
% this year had strong La Nina event
figure;
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(vard(:,ii+12*9),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
caxis([-30 30])
hold on
plot(long,lat,'k');
plot(long+360.,lat,'k');
ylim([45 60]); 
xlim([220 245])
if ii==12
colorbar
end
title([months{ii}, ' 1988'])
end

% plot temp for all months of 1997
% this year had strong El Nino event
figure;
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(vard(:,ii+12*18),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
caxis([-30 30])
hold on
plot(long,lat,'k');
plot(long+360.,lat,'k');
ylim([45 60]); 
xlim([220 245])
if ii==12
colorbar
end
title([months{ii}, ' 1997'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing: 
% remove the seasonal cycle and apply 3-month running mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% removing 7 months in 2017 (so the timeseries is for 1979-2017)
vard_cut=vard(:,1:end-7);

nyear=38;
% calculate seasonal cycle for each grid point
for i=1:length(vard(:,1))
var_seasonal(i,:)=nanmean((reshape(vard_cut(i,:),12,nyear))');
end

% copy-paste this cycle for all the years
var_seasonal_all=repmat(var_seasonal,1,nyear);
% add Jan-Jul for 2017 so the seasonal cycle timeseries is the same lenght as original timeseries
var_seasonal_all=[var_seasonal_all var_seasonal(:,1:7)];
% remove the seasonal cycle from the original series 
vard_new=vard-var_seasonal_all;

% apply 3-month running mean
window=3;
for i=1:length(vard(1,:))-window+1;
    runmean(i,:)=nanmean(vard_new(:,i:i+window-1)');
end
vard_runmean=runmean';

% lets plot all these pre-processing steps for one grid point (grid point number 50):
% timeseries for original data
timen=[1979+1/12:1/12:2017+7/12];
% timeseries for the running mean data (since it has 2 fewer points than the original series
% start with Feb 1997 and end with May 2017
time=[1979+2/12:1/12:2017+6/12];

figure;
subplot(4,1,1)
plot(timen,vard(50,:),'b-');
xlabel('Time (yr)');
ylabel('T (^oC)');
title('original time series for one grid point')
xlim([1979 2017.6]);

subplot(4,1,2)
plot(timen,var_seasonal_all(50,:),'b-');
xlabel('Time (yr)');
ylabel('T (^oC)');
title('average seasonal for the grid point')
xlim([1979 2017.6]);

subplot(4,1,3)
plot(timen,vard_new(50,:),'b-');
xlabel('Time (yr)');
ylabel('T anomaly (^oC)');
title('residual (original-seasonal cycle)')
xlim([1979 2017.6]);

subplot(4,1,4)
plot(time,vard_runmean(50,:),'b-');
xlabel('Time (yr)');
ylabel('T anomaly (^oC)');
title('3-month running mean of the residual')
xlim([1979 2017.6]);

% plot the 3-month running mean of residuals for all
% grid point for all months of 1988
figure;
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(vard_runmean(:,ii+12*9-1),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
caxis([-3 3])
hold on
plot(long,lat,'k');
plot(long+360.,lat,'k');
ylim([45 60]); 
xlim([220 245])
if ii==12
colorbar
end
title(['runmean ', months{ii}, ' 1988'])
end

% plot the 3-month running mean of residuals for all
% grid point for all months of 1997
figure;
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(vard_runmean(:,ii+12*18-1),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
caxis([-3 3])
hold on
plot(long,lat,'k');
plot(long+360.,lat,'k');
ylim([45 60]); 
xlim([220 245])
if ii==12
colorbar
end
title(['runmean ',months{ii}, ' 1997'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xdata is SST anomalies (Feb 1979 to May 2017)
% ydata will be SLP anomalies (need to cut one month since it goes until Jun 2017)
xdata=xdata;
data=vard_runmean';  % BC temp data is the 3-mont running mean calculated above
ydata=data(1:end-1,:);

[A,B,r,U,V] = canoncorr(xdata,ydata);

% to map the CCA modes back to the original space:
% x=FU, y=GV
% where F=cov(xdata)*A; and G=cov(ydata)*B;

F=cov(xdata)*A;
G=cov(ydata)*B;

time=[1979+2/12:1/12:2017+5/12];

% plot first three CCA modes (note that 450 modes are produces) 
br=0;
figure; 
for i=1:3
br=br+1;
subplot(3,3,br)
var01=F(:,i);
var_plot=reshape(var01,tt_SST(1),tt_SST(2));
h=imagesc(x_SST,y_SST,var_plot');
axis xy
colormap jet
colorbar;
title(['SST CCA mode ',num2str(i)]);

br=br+1;
subplot(3,3,br)
var01=G(:,i);
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
colorbar;
hold on
plot(long,lat,'k');
plot(long+360.,lat,'k');
ylim([45 60]); 
xlim([220 245])
title(['temp CCA mode ',num2str(i)]);

br=br+1;
subplot(3,3,br)
plot(time,U(:,i),'r-',time,V(:,i),'b-')
xlabel('time');
xlim([1979 2017.6]);
title(['r=',num2str(r(i),'%2.2f')]);
if i==1
legend('u(SST)','v(temp)');
end
end


% reconstructing all xdata from the first PCA modes
% I decided to keep the first 2 modes for SST data
[xeigenvectors,xPCs,xeigenvalues]=princomp(xdata);
invE=inv(xeigenvectors');
xPCs=xPCs';
xdata_rec=(invE(:,1:2)*xPCs(1:2,:))';

% reconstructing all xdata from the first PCA modes
% I decided to keep the first 5 modes for BC temp data
[yeigenvectors,yPCs,yeigenvalues]=princomp(ydata);
invE=inv(yeigenvectors');
yPCs=yPCs';
ydata_rec=(invE(:,1:5)*yPCs(1:5,:))';


% CCA applied on pre-processed data 
[A,B,r,U,V] = canoncorr(xdata_rec,ydata_rec);

F=cov(xdata)*A;
G=cov(ydata)*B;

time=[1979+2/12:1/12:2017+5/12];

% plot first two CCA modes (note that only two are produced)  
br=0;
figure; 
for i=1:2
br=br+1;
subplot(2,3,br)
var01=F(:,i);
var_plot=reshape(var01,tt_SST(1),tt_SST(2));
h=imagesc(x_SST,y_SST,var_plot');
axis xy
colormap jet
colorbar;
title(['Reconstr.SST CCA mode ',num2str(i)]);

br=br+1;
subplot(2,3,br)
var01=G(:,i);
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
colorbar;
hold on
plot(long,lat,'k');
plot(long+360.,lat,'k');
ylim([45 60]); 
xlim([220 245])
title(['Reconstr.temp CCA mode ',num2str(i)]);

br=br+1;
subplot(2,3,br)
plot(time,U(:,i),'r-',time,V(:,i),'b-')
xlabel('time');
xlim([1979 2017.6]);
title(['r=',num2str(r(i),'%2.2f')]);
if i==1
legend('Reconstr. u(SST)','Reconstr. v(temp)');
end
end
