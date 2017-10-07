%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 4 (26 Sep 2017)
% PCA on real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2:
% Gridded monthly sea surface temperature data for Tropical Pacific
% from ERA Interim reanalysis
% Period: Jan 1979 to Jun 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

months={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};

% this is nc formated file (binary, structured file)
filename='SST_ERAInterim_monthly_Jan1979_Jun2017_Tropical_Pacific.nc'
% the following command will display the file content
ncdisp(filename);

% these are the factors read from the file (should be displayed by 'ncdisp');
scale_factor= 0.00018735;
add_offset= 298.0465;
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
C_var=strcmp(varname,'sst');
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
% var will have sst values
var = netcdf.getVar(ncid,varid);
end

end

% closing the nc file
netcdf.close(ncid)

tt=size(var);
% set missing values to NaN
var(var == -32767) = NaN;
% set sst to deg C (and scale it accoring to the parameters given in the file)
var=var.*scale_factor + add_offset -273.15;  

% var is 3-D matrix, here I am reshaping it into 2-D
% where rows are grid points, columns are time points
for kj=1:tt(3)
var1=squeeze(var(:,:,kj));
vard(1:tt(1)*tt(2),kj)=reshape(var1,tt(1)*tt(2),1);
end

% converting integer values to double precision 
vard=double(vard);

% lets plot the sst over the whole domain for Jan 1979
var_plot=reshape(vard(:,1),tt(1),tt(2));
%load coast 
% if you don't have the map toolbox use coastline.mat:
load 'coastline.mat';

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

% plot sst for all months in 1988
% this year had strong La Nina event
figure;
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(vard(:,ii+12*9),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
caxis([22 31])
if ii==12
colorbar
end
title([months{ii}, ' 1988'])
end

% plot sst for all months of 1997
% this year had strong El Nino event
figure;
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(vard(:,ii+12*18),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
caxis([22 31])
if ii==12
colorbar
end
title([months{ii}, ' 1997'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing: 
% remove seasonal cycle and apply 3-month running mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% removing the 6 months in 2017 (so the timeseries is for Jan 1979 to Dec 2016)
% total of 38 years
nyear=38;

vard_cut=vard(:,1:end-6);

% calculate a seasonal cycle for each grid point
for i=1:length(vard(:,1))
var_seasonal(i,:)=nanmean((reshape(vard_cut(i,:),12,nyear))');
end

% copy-paste this cycle for all the years
var_seasonal_all=repmat(var_seasonal,1,nyear);
% add Jan-Jun for 2017 so the seasonal cycle timeseries is the same lenght as original timeseries
var_seasonal_all=[var_seasonal_all var_seasonal(:,1:6)];
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
timen=[1979+1/12:1/12:2017+6/12];
% timeseries for the running mean data (since it has 2 fewer points than the original series
% start with Feb 1997 and end with May 2017
time=[1979+2/12:1/12:2017+5/12];

figure;
subplot(4,1,1)
plot(timen,vard(50,:),'b-');
xlabel('Time (yr)');
ylabel('SST (^oC)');
title('original time series for one grid point')
xlim([1979 2017.6]);
ylim([24 30]);

subplot(4,1,2)
plot(timen,var_seasonal_all(50,:),'b-');
xlabel('Time (yr)');
ylabel('SST (^oC)');
title('average seasonal for the grid point')
xlim([1979 2017.6]);
ylim([24 30])

subplot(4,1,3)
plot(timen,vard_new(50,:),'b-');
xlabel('Time (yr)');
ylabel('SST anomaly (^oC)');
title('residual (original-seasonal cycle)')
xlim([1979 2017.6]);

subplot(4,1,4)
plot(time,vard_runmean(50,:),'b-');
xlabel('Time (yr)');
ylabel('SST anomaly (^oC)');
title('3-month running mean of the residual')
xlim([1979 2017.6]);

% plot the 3-month running mean of residuals for all
% grid points for all months of 1988
figure;
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(vard_runmean(:,ii+12*9-1),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
caxis([-3 3])
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
if ii==12
colorbar
end
title(['runmean ',months{ii}, ' 1997'])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preparing the data for PCA analysis
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
xlim([1979 2017.6]);
end

% finding El Nino and La Nina events
[pks1,indEN]=findpeaks(PCs(:,1),'MinPeakHeight',45,'MinPeakWidth',0);
[pks2,indLN]=findpeaks(-PCs(:,1),'MinPeakHeight',45,'MinPeakWidth',2);

indEN1=indEN(pks1 >=80); 
indLN1=indLN(pks2 >=85);

figure; 
plot(time,PCs(:,1),'k-'); hold on
plot(time(indEN),PCs(indEN,1),'ro'); hold on
plot(time(indLN),PCs(indLN,1),'bo'); hold on
plot(time(indEN1),PCs(indEN1,1),'r*'); hold on
plot(time(indLN1),PCs(indLN1,1),'b*'); hold on
xlabel('time');
title(['PC',num2str(i)]);
xlim([1950 2017.6]);
grid on

% plot data in space of the two first two PCs
figure;
plot(PCs(:,1),PCs(:,2),'k.',PCs(indEN,1),PCs(indEN,2),'r.',PCs(indLN,1),PCs(indLN,2),'b.'); hold on
plot(PCs(indEN1,1),PCs(indEN1,2),'r*',PCs(indLN1,1),PCs(indLN1,2),'b*');
xlabel('PC1');
ylabel('PC2');
legend('all','strong El Nino','strong La Nina');
grid on



