%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 9 (31 Oct 2017)
% Application of SOM on a synthetic data 
% Example 2 from the tutorial

% NOTE: in order to run this script you need to have
% SOM Toolbox (downloaded from the website:
% http://www.cis.hut.fi/projects/somtoolbox/)
% included in your Matlab path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creating the input data
k=2*pi/100;
omega=2*pi/50;
x=[1:100];
t=[1:200];

for i=1:200
y(i,:)=sin(k*x-omega*t(i));
end

% sine pattern w/ amplitude of 1
y1=y(1,:);

% cosine pattern w/ amplitude of 0.5
y2=0.5*y(37,:);

% step pattern w/ amplitude of 0.8
y3(1:50)=-0.8;
y3(51:100)=0.8;

% sawtooth pattern w/ amplitude of 1
y4(1:25)=-2*x(1:25)./25+1;
y4(26:50)=2*x(26:50)./25-3;
y4(51:75)=-2*x(51:75)./25+5;
y4(76:100)=2*x(76:100)./25-7;

% create timeseries where from t=1:50 signal y1 occurs,
% then for the next 50 steps signal y3 occurs,
% then for the next 50 steps signal y4 occurs,
% and then for the last 50 steps signal y2 occurs
ynew(1:50,:)=repmat(y1,50,1);
ynew(51:100,:)=repmat(y3,50,1);
ynew(101:150,:)=repmat(y4,50,1);
ynew(151:200,:)=repmat(y2,50,1);

data=ynew;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOM algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initilizing SOM
% chose size of the map for SOM 
% I tried 2x2, 2x3 and 3x2
ny_som=2; nx_som=2;
en=ny_som*nx_som;

data=double(data);

msize=[ny_som nx_som];
% performing linear initialization of nodes
display('initialization')
sMap=som_lininit(data,'msize',msize,'hexa','sheet');

% training SOM
display('training')
[sM,sT] = som_batchtrain(sMap,data,'ep','hexa','sheet','radius',[2 1],'trainlen',200); 

% calulating quantization and topological error
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
title(['SOM node' num2str(index(i)) ', f=' num2str(hi(index(i)),'%2.1f')])
end

% derive timeseries of SOM patterns
bmus=som_bmus(sM,data);
figure; 
plot(bmus);
xlabel('time');
ylabel('SOM node index');

% some additional plots:
%% plot colored SOM (with actual hexagonal nodes)
%  and distance among neighboring nodes 
U = som_umat(sM);
Um = U(1:2:size(U,1),1:2:size(U,2));
C = som_colorcode(sM,'rgb2');

figure('Renderer','Painters')
som_cplane(sM,C,1-Um(:)/max(Um(:)));
title('Color coding + distance matrix')

% plot colored SOM (with actual hexagonal nodes)
% and visualize the SOM patterns inside each node
figure('Renderer','Painters')
som_cplane(sM,C);
hold on
som_plotplane(sM,sM.codebook)
title('SOM with patterns')
set(gca,'xticklabel',[])
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDITIONAL material 
% the following algorithm determines optimal SOM size, however it's usualy a very large SOM
sD = som_data_struct(data); 
tic
sM = som_make(sD); 
toc
ny_som = sM.topol.msize(1);
nx_som = sM.topol.msize(2);
en=ny_som*nx_som;

bmus=som_bmus(sM,data);

figure; 
plot(bmus);
xlabel('time');
ylabel('SOM node index');

myhits = som_hits(sM,sD);
[inbNaN dummy]=find(myhits == 0);


%% plot colored SOM (with actual hexagonal nodes)
%  and distance among neigbouring nodes 
U = som_umat(sM);
Um = U(1:2:size(U,1),1:2:size(U,2));
C = som_colorcode(sM,'rgb2');

C0=C;
C0(inbNaN,1)=1;
C0(inbNaN,2)=1;
C0(inbNaN,3)=1;

figure('Renderer','Painters')
som_cplane(sM,C,1-Um(:)/max(Um(:)));
title('Optimized SOM with distances')

% plot colored SOM (with actual hexagonal nodes)
% and SOM patterns
figure('Renderer','Painters')
som_cplane(sM,C0);
hold on
som_plotplane(sM,sM.codebook)
title('Optimized SOM with patterns')
set(gca,'xticklabel',[])
hold off


