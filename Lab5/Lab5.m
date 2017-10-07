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
tt_SST=[135 29 462];
x_SST=[159.75:0.75:260.25]';
y_SST=[10.50:-0.75:-10.50]';
% put this data into variable xdata
xdata=data_SST;
clear data

months={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};

% this is nc formated file (binary, structured file)
filename='T_monthly_Canada_1979_2017_ERA_Interim.nc';

% the following command will display the file content
ncdisp(filename);

% these are the factors read from the file
scale_factor  = 0.001107
add_offset    = 271.712
missing_value = -32767
%units         = 'K'
%long_name     = '2 metre temperature'

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
C_var=strcmp(varname,'t2m');
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
% remove seasonal cycle and apply 3-month running mean
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
% start with Feb 1997 and end with May 2015
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


