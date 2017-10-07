%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 5 (3 Oct 2017)
% CCA on the synthetic data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load some data
filename='Tutorial5_input_data.mat';
load(filename);

[n1 n2]=size(data);
% create the input data (xdata and ydata)
xdata=data(:,1:3);
ydata=2*data(:,[4 1]);

% time variable:
t=[1:n1];

% run the code first with no rotation in ydata 
% (set theta to zero)
% then run with the rotation of 30 deg counterclockwise
theta=-30;
R=[cosd(theta) sind(theta); -sind(theta) cosd(theta)];

ydata=R*ydata';

ydata=ydata';

% plot the timeseries for each variable
figure;
subplot(3,2,1);
plot(t,xdata(:,1),'b.-','Linewidth',1);
ylabel('x1');
title('x data');

subplot(3,2,3);
plot(xdata(:,2),'b.-','Linewidth',1);
ylabel('x2');

subplot(3,2,5);
plot(xdata(:,3),'b.-','Linewidth',1);
ylabel('x3');
xlabel('time');

subplot(3,2,2);
plot(ydata(:,1),'r.-','Linewidth',1);
ylabel('y1');
title('y data')

subplot(3,2,4);
plot(ydata(:,2),'r.-','Linewidth',1);
ylabel('y2');
xlabel('time');

% these are variables for scatter plot only
S(1:n1)=10;
C=[1:n1];
C=C';

% scatter plots
figure;
subplot(1,2,1);
scatter3(xdata(:,1),xdata(:,2),xdata(:,3),S,C,'filled');
xlim([-400 200]);
ylim([-100 50]);
zlim([-30 40]);
xlabel('x1');
ylabel('x2');
zlabel('x3');
colorbar('SouthOutside');
colormap jet 
title('x data')

subplot(1,2,2);
scatter(ydata(:,1),ydata(:,2),S,C,'filled');
xlim([-80 60]);
ylim([-100 100]);
xlabel('y1');
ylabel('y2');
colorbar('SouthOutside');
colormap jet
title('y data');

% CCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A,B,r,U,V] = canoncorr(xdata,ydata);

% to map the CCA modes back to the original space:
% x=FU, y=GV
% where F=cov(X)*A and G=cov(Y)*B

F=cov(xdata)*A;
G=cov(ydata)*B;

% plot the 2 CCA modes (vectors F and G, and temporal coefficients U and V)
br=0;
figure; 
for i=1:2
br=br+1;
subplot(2,3,br)
plot([1:3],F(:,i),'bo-');
if i==1
title('x mode')
end

br=br+1;
subplot(2,3,br)
plot([1:2],G(:,i),'ro-');
if i==1
title('y mode')
end

br=br+1;
subplot(2,3,br)
plot([1:n1],U(:,i),'b-',[1:n1],V(:,i),'r-');
xlabel('time');
if i==1
legend('u','v');
end
end

S(1:n1)=10;
C=[1:n1];
C=C';

% scatter plots with vectors F and G
figure;
subplot(1,2,1);
scatter3(xdata(:,1),xdata(:,2),xdata(:,3),S,C,'filled'); hold on
plotv([8*F(:,1) 8*F(:,2)],'-');
xlim([-400 200]);
ylim([-100 50]);
zlim([-30 40]);
xlabel('x1');
ylabel('x2');
zlabel('x3');
legend('data','F1','F2');
colorbar('SouthOutside');
colormap jet 
title('x data')

subplot(1,2,2);
scatter(ydata(:,1),ydata(:,2),S,C,'filled'); hold on
plotv([3*G(:,1) 5*G(:,2)],'-');
xlim([-80 60]);
ylim([-100 100]);
xlabel('y1');
ylabel('y2');
legend('data','G1','G2');
colorbar('SouthOutside');
colormap jet
title('y data');


