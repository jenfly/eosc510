%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 9 (2 Nov 2016)
% Application of SOM on a synthetic data 
% Example 3 from the tutorial
% linear progressive wave y(x,t)=sin(k*x-omega*t)
% NOTE: in order to run this script you need to have
% SOM Toolbox installed from the website:
% http://www.cis.hut.fi/projects/somtoolbox/
% in your Maltab path include the downloaded folder with som toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creating the input data
k=2*pi/100;
omega=2*pi/50;
x=[1:100];
t=[1:200];

for i=1:200
y(i,:)=sin(k*x-omega*t(i));
end

% plot spatial wave (in x domain) for each time step (t=1,...,200)
figure;
for i=1:50
subplot(5,10,i)
plot(y(i,:));
ylim([-1 1])
title(['t=',num2str(i)])
end

% plot temporal wave (in t domain) for each point in x (x=1,...100)
figure;
for i=1:50
subplot(5,10,i)
plot(y(:,i));
ylim([-1 1])
title(['x=',num2str(i)])
end

data=y;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOM algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initilizing SOM
% chose size of the map for SOM 
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

% plot the inital SOM patterns
figure;
for i=1:en
subplot(ny_som,nx_som,i);
plot([1:100],sMap.codebook(index(i),:),'b-');
xlabel('time');
title(['Initial node' num2str(index(i))])
end

% plot the final (after training) SOM patterns
figure;
for i=1:en
subplot(ny_som,nx_som,i);
plot([1:100],sM.codebook(index(i),:),'b-');
xlabel('time');
title([num2str(index(i)) ', f=' num2str(hi(index(i)),'%2.1f')])
end

% derive timeseries of SOM patterns
bmus=som_bmus(sM,data);
figure; 
plot(bmus);
xlabel('time');
ylabel('SOM node index');
