%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 5 (3 Oct 2017)
% CCA on real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we take the SST anomalies data from last tutorial (saved in SST_anomalies_data.mat):
% and add the same resolution gridded data of
% monthly mean sea level pressure over the same domain (Tropical Pacific)
% Period: Jan 1979 to July 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load SST anomalies (the processed data from last Tutorial)
load SST_anomalies_data.mat;
% put this data into variable xdata
xdata=data;
clear data

months={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};

% this is nc formated file (binary, structured file)
filename='SLP_ERAInterim_monthly_Jan1979_Jul2017_Tropical_Pacific.nc';

% the following command will display the file content
ncdisp(filename);

% these are the factors read from the file
scale_factor  = 0.016543
add_offset    = 101013.8433
missing_value = -32767;

% opening the file
ncid = netcdf.open(filename,'NC_NOWRITE');   % this is the matlab function for opening nc files
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);  % this is all the structure in the nc file

% this loop looks for a particular variable from the set of variables 
for jj=1:numvars
[varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,jj-1);
varname=varname % display name of the variables

% strcmp is a logical function used here for finding each on for the following variables: lon, lat, sst, time
% we know the names from varname and now want to pull out the values for these variables
C_lon=strcmp(varname,'longitude'); 
C_lat=strcmp(varname,'latitude'); 
C_var=strcmp(varname,'msl');
C_time=strcmp(varname,'time');

if C_lon == 1
% Get variable ID of the first variable, given its name.
varid = netcdf.inqVarID(ncid,varname);
% assign that variable values to x
% so x will be longitudes (in degrees East)
x = netcdf.getVar(ncid,varid);
end

if C_lat == 1
varid = netcdf.inqVarID(ncid,varname);
% y will be latitudes in degrees North
y = netcdf.getVar(ncid,varid);
end

if C_time == 1 
varid = netcdf.inqVarID(ncid,varname);
time = netcdf.getVar(ncid,varid);
end

if C_var == 1 
varid = netcdf.inqVarID(ncid,varname);
% var will have msl (mean sea level) values
var = netcdf.getVar(ncid,varid);
end

end

% closing the nc file
netcdf.close(ncid)

tt=size(var);

% set missing values to NaN
var(var == -32767) = NaN;
% set sst to deg C (and scale it accoring to the parameters given in the file)
var=double(var);
var=var.*scale_factor + add_offset;
var=var./100; % setting units to hPa  

% var is 3-D matrix, here I am reshaping it into 2-D
% where rows are grid points, columns are time points
for kj=1:tt(3)
var1=squeeze(var(:,:,kj));
vard(1:tt(1)*tt(2),kj)=reshape(var1,tt(1)*tt(2),1);
end

% converting integer values to double precision 
vard=double(vard);

% lets plot the msl over the whole domain for Jan 1979
var_plot=reshape(vard(:,1),tt(1),tt(2));
load coast
figure;
h=imagesc(x,y,var_plot');
colormap jet
colorbar
axis xy
hold on
plot(long,lat,'k');
plot(long+360.,lat,'k');
ylim([-50 50]); 
xlim([120 300])
hold off

% plot msl for all months in 1988
% this year had strong La Nina event
figure;
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(vard(:,ii+12*9),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
caxis([1008 1014])
if ii==12
colorbar
end
title([months{ii}, ' 1988'])
end

% plot msl for all months of 1997
% this year had strong El Nino event
figure;
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(vard(:,ii+12*18),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
caxis([1008 1014])
if ii==12
colorbar
end
title([months{ii}, ' 1997'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing: 
% remove seasonal cycle and apply 3-month running mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% removing 7 months in 2017 (so the timeseries is for 1979-2014)
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
% start with Feb 1997 and end with May 2015
time=[1979+2/12:1/12:2017+6/12];

figure;
subplot(4,1,1)
plot(timen,vard(50,:),'b-');
xlabel('Time (yr)');
ylabel('SLP (hPa)');
title('original time series for one grid point')
xlim([1979 2017.6]);
ylim([1005 1015]);

subplot(4,1,2)
plot(timen,var_seasonal_all(50,:),'b-');
xlabel('Time (yr)');
ylabel('SLP (hPa)');
title('average seasonal for the grid point')
xlim([1979 2017.6]);
ylim([1005 1015]);

subplot(4,1,3)
plot(timen,vard_new(50,:),'b-');
xlabel('Time (yr)');
ylabel('SLP anomaly (hPa)');
title('residual (original-seasonal cycle)')
xlim([1979 2017.6]);

subplot(4,1,4)
plot(time,vard_runmean(50,:),'b-');
xlabel('Time (yr)');
ylabel('SLP anomaly (hPa)');
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
caxis([-2 2])
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
caxis([-2 2])
if ii==12
colorbar
end
title(['runmean ',months{ii}, ' 1997'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparing the data for PCA analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=vard_runmean';

% PCA:
[eigenvectors,PCs,eigenvalues]=princomp(data);

% contribution of each mode to total variance:
variance=eigenvalues./sum(eigenvalues);

figure; 
subplot(2,1,1)
plot([1:length(variance)],variance,'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

subplot(2,1,2)
plot([1:20],variance(1:20),'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

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
xlim([1979 2017.7]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xdata is SST anomalies (Feb 1979 to May 2017)
% ydata is SLP anomalies (cut one month since it goes until Jun 2017)
xdata=xdata;
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
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
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
title(['SLP CCA mode ',num2str(i)]);

br=br+1;
subplot(3,3,br)
plot(time,U(:,i),'r-',time,V(:,i),'b-')
xlabel('time');
xlim([1979 2017.6]);
title(['r=',num2str(r(i),'%2.2f')]);
if i==1
legend('u(SST)','v(SLP)');
end
end

% reconstruct the data using only first 2 modes for both dataset
% reconstructing all ydata
for j=1:length(ydata(:,1))
ydata_rec(j,:)=PCs(j,1).*eigenvectors(:,1)+PCs(j,2).*eigenvectors(:,2);
end

[xeigenvectors,xPCs,xeigenvalues]=princomp(xdata);
% reconstructing all xdata
for j=1:length(xdata(:,1))
xdata_rec(j,:)=xPCs(j,1).*xeigenvectors(:,1)+xPCs(j,2).*xeigenvectors(:,2);
end

% plot the 3-month running mean of residuals for all
% grid point for all months of 1997
figure;
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(xdata_rec(ii+12*18-1,:),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
caxis([-3 3])
if ii==12
colorbar
end
title(['SST rec ',months{ii}, ' 1997'])
end

figure;
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(ydata_rec(ii+12*18-1,:),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
caxis([-2 2])
if ii==12
colorbar
end
title(['SLP rec ',months{ii}, ' 1997'])
end

% CCA applied on pre-processed data (data reconstructed only from first 3 modes)
% xdata is SST anomalies (Feb 1979 to May 2017)
% ydata is SLP anomalies (cut one month since it goes until Jun 2017)
[A,B,r,U,V] = canoncorr(xdata_rec,ydata_rec);

% to map the CCA modes back to the original space:
% x=FU, y=GV
% where F=cov(X)*A and G=cov(Y)*B

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
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
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
title(['Reconstr.SLP CCA mode ',num2str(i)]);

br=br+1;
subplot(2,3,br)
plot(time,U(:,i),'r-',time,V(:,i),'b-')
xlabel('time');
xlim([1979 2017.6]);
title(['r=',num2str(r(i),'%2.2f')]);
if i==1
legend('Reconstr. u(SST)','Reconstr. v(SLP)');
end
end

