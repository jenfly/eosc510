%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 9 (2 Nov 2016)
% SOMs on real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2:
% This is the same data as from Tutorial 4 on PCA:
% Gridded monthly sea surface temperature data for Tropical Pacific
% from ERA Interim reanalysis
% Period: Jan 1979 to Jun 2015
% The first part of the script is a copy from Tuturial4_SLP_data.m
% NOTE: in order to run this script you need to have
% SOM Toolbox installed from the website:
% http://www.cis.hut.fi/projects/somtoolbox/
% in your Maltab path include the downloaded folder with som toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

months={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};

% this is nc formated file (binary, structured file)
% create your path here for the directory where the data from Tutorial 4 is (these paths are for my computer):
filename='/shared/EOSC510/Tutorials/Tutorial4/Example2/SST_monthly_1979_2015_ERAInterim.nc';

% the following command will display the file content
ncdisp(filename);

% these are the factors read from the file
scale_factor= 0.00018739;
add_offset= 298.0479;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing: 
% remove seasonal cycle and apply 3-month running mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% removing 6 months in 2015 (so the timeseries is for 1979-2014)
vard_cut=vard(:,1:end-6);

% calculate seasonal cycle for each grid point
for i=1:length(vard(:,1))
var_seasonal(i,:)=nanmean((reshape(vard_cut(i,:),12,36))');
end

% copy-paste this cycle for all 36 years
var_seasonal_all=repmat(var_seasonal,1,36);
% add Jan-Jun for 2015 so the seasonal cycle timeseries is the same lenght as original timeseries
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
timen=[1979+1/12:1/12:2015+6/12];
% timeseries for the running mean data (since it has 2 fewer points than the original series
% start with Feb 1997 and end with May 2015
time=[1979+2/12:1/12:2015+5/12];

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

% extend bmus with NaN so that is starts with Jan 1979 to Dec 2015
bmusn(1:37*12)=NaN;
bmusn(2:end-7)=bmus;

% plot bmus in 2D
x1=[1:12];
y1=[1979:2015];
figure;
bmus_plot=reshape(bmusn,12,37);
h=imagesc(y1,x1,bmus_plot);
colormap jet
colorbar
set(h,'alphadata',~isnan(bmus_plot));
xlabel('Time (yr)');
ylabel('Month');


