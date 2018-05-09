%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 9 (31 Oct 2017)
% Application of SOM on a synthetic data 
% Example 3 from the tutorial
% linear progressive wave y(x,t)=sin(k*x-omega*t)

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




