%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 2 (14 Sep 2016)
% loading Crime_data.csv
% perform multiple linear regression and stepwise regression
% rank the importance of the individual predictors in their influence on y
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data=importdata('Crime_data.csv');

% Column X1 and X2 are the response variables 
% y, first column: total overall reported crime rate per 1 million residents
% y, second column: reported violent crime rate per 100,000 residents

y=data.data(:,2); % setting y to be equal to 1st or 2nd column

% other columns (X3 to X7) are predictors
xall=data.data(:,[3:end]); % setting the rest of the columns to be matrix X

n=size(y,1);

% MLR with original variables
X=[ones(n,1) xall];   % first column needs to consist of ones
b=regress(y,X)
y_regr=X*b;

% stepwise regression with original values
[a SE PVAL INMODEL STATS NEXTSTEP HISTORY]=stepwisefit(xall,y,'penter',0.05);

% regression coefficients
a=a

% constant coefficient 
a0=STATS.intercept

% in order to rank the importance of individual predictors lets first check whether 
% we need to standardize them (i.e. how much their mean and std differ)
mean_xall=mean(xall)
std_xall=std(xall)

% standardizing the predictors (x1-mean(x1))/std(x1), etc
xall_mean=repmat(mean_xall,n,1);
xall_std=repmat(std_xall,n,1);

xall_standard=(xall-xall_mean)./xall_std;

%regression with standardized variables
X=[ones(n,1) xall_standard];   
[b2,BINT2,R2,RINT2,STATS2]=regress(y,X);
y_regr=X*b2;

% stepwise regression with standardized values
[a2 SE PVAL INMODEL STATS NEXTSTEP HISTORY]=stepwisefit(xall_standard,y,'penter',0.05);

% regression coefficients
a2=a2

% constant coefficient 
a20=STATS.intercept

% plot the standardized variables (and standardized y)
y_standard=(y-mean(y))/std(y);

figure;
subplot(2,1,1);
plot([1:50],y,'LineWidth',1);
legend('# crimes');

subplot(2,1,2);
plot([1:50],xall,'LineWidth',1);
legend('police funding','25+ w/ 4 yrs of highschool','16-19 not in highschool','18-25 in college','25+ w/ 4 yrs in collegue');

figure;
plot([1:50],[y_standard xall_standard],'LineWidth',1);
legend('# crimes','police funding','25+ w/ 4 yrs of highschool','16-19 not in highschool','18-25 in college','25+ w/ 4 yrs in collegue');
title('X standardized')

% take out only the relevant predictors (from stepwisefit)
[ind1 ind2]=find(INMODEL == 1);
b=regress(y,X(:,[1 ind2+1]))
y_regr=X(:,[1 ind2+1])*b;

% plot Y_regr vs Y and provide correlation between them
[r p]=corrcoef([y y_regr]);

xline=[min([y;y_regr]):max([y;y_regr])]; 
figure;
plot(y,y_regr,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([y;y_regr]) max([y;y_regr])]);
ylim([min([y;y_regr]) max([y;y_regr])]);
xlabel('y');
ylabel('y_{regr}')
title(['Stepwise regression r= ' num2str(r(1,2),'%2.2f')]);

% check individual correlations (y and each column of X):
[r p]=corrcoef([y xall])
