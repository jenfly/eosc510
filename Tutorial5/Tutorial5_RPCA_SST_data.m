%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 5 (3 Oct 2017)
% Rotated PCA on real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example on SST data (from the last Tutorial)
% Gridded monthly sea surface temperature data for Tropical Pacific
% from ERA Interim reanalysis
% Period: Jan 1979 to Jun 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load SST anomalies (the processed data from last Tutorial)
load SST_anomalies_data.mat;
% loads variable 'data'

% a variable needed for plots
tt=[135 29 462];
x=[159.75:0.75:260.25]';
y=[10.50:-0.75:-10.50]';

% timeseries for the running mean data 
% start with Feb 1997 and end with May 2017
time=[1979+2/12:1/12:2017+5/12];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotated PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% these variable will be needed for plotting eigenvectors (coordinates)
X1=[0 max(PCs(:,1))];
Y1=[0 max(PCs(:,2))];
Z1=[0 max(PCs(:,3))];

% here we apply the rotation (varimax criterion) for the first three PCs
% decreased the number of time points to 450 (460 was too close to the iteration limit for the rotation)
[RPCs,R] = rotatefactors(PCs(1:450,1:3));

% plot data in space of the first two modes
% plot eigenvectors (coordinates) with thin lines (Linewidth:1)
% plot rotated eignevectors with thick lines (Linewidth: 2)

% coordinates for rotated values are given in matrix R  
% So for rotated eigenvector 1: R(1,1)i+R(2,1)j  (where i and j are unit vectors)
% For rotated eigenvector 2: R(1,2)i+R(2,2)j

figure;
plot(PCs(:,1),PCs(:,2),'k.',PCs(indEN,1),PCs(indEN,2),'r.',PCs(indLN,1),PCs(indLN,2),'b.'); hold on
plot(PCs(indEN1,1),PCs(indEN1,2),'r*',PCs(indLN1,1),PCs(indLN1,2),'b*'); hold on
line(X1,[0 0],'Color','r','LineWidth',1);
line([0 0],Y1,'Color','g','LineWidth',1);
line([0 R(1,1)*100],[0 R(2,1)*100],'Color','r','LineWidth',2);
line([0 R(1,2)*100],[0 R(2,2)*100],'Color','g','LineWidth',2);
xlabel('PC1');
ylabel('PC2');
legend('all','El Nino','La Nina');
grid on
xlim([-100 150])
ylim([-30 60])

% find rotated eigenvectors
Reigenvectors=eigenvectors(:,1:3)*R;

% calulate how much variance each RPCA mode carries
% total variance
var_total=sum(diag(cov(data)));
% lambdas (eigenvalues) of RPCA
for i=1:3
lambda(i)=Reigenvectors(:,i)'*cov(data)*Reigenvectors(:,i);
end
% explained variance for RPCA
explained=lambda./var_total;

% plot first 3 modes for eignevectors, rotated eigenvectors, PCs and RPCs
br=0;
figure; 
for i=1:3
br=br+1;
subplot(3,3,br)
var01=eigenvectors(:,i);
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
colorbar;
title(['eigenvector var=',num2str(variance(i),'%2.2f')]);

br=br+1;
subplot(3,3,br)
var01=Reigenvectors(:,i);
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
colorbar;
title(['Reigenvector var=',num2str(explained(i),'%2.2f')]);

br=br+1;
subplot(3,3,br)
plot(time(1:450),PCs(1:450,i),'b-',time(1:450),RPCs(1:450,i),'r-')
xlabel('time');
title(['PC',num2str(i)]);
xlim([1979 2017.6]);
if i==1
legend('PC','RPC');
end
end


