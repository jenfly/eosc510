%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 2 (14 Sep 2016)
% script that goes with the presentation from the Tutorial
% Four examples on applying multiple linear regression and stepwise regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1: testing MLR on artificial data
% Y=a0+a1*X1+a2*X2+a3*X3+a4*X4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load X data (created as X=randi(100,50,4))
% X=[X1 X2 X3 X4];

load Xdata.mat;

% define the regression coefficients
a0=0;
a1=1;
a2=-2;
a3=3;
a4=-4;

Y=a0+a1*X(:,1)+a2*X(:,2)+a3*X(:,3)+a4*X(:,4);

% plot Y and X
n=size(Y,1);
figure;
plot([1:n],Y,'k-','LineWidth',2);
hold on
plot([1:n],X,'LineWidth',1);
legend('Y','X1','X2','X3','X4');

% applying MLR on Y and X
n=size(Y,1);
Xnew=[ones(n,1) X];   % first column needs to consist of ones
b=regress(Y,Xnew)
Y_regr=Xnew*b;

% Instead of using regress.m lets derive the coefficients from 
% minimizing sum of squared errors (SSE) with respect to a 
%a=(X^T X)^-1 X^T y

areg=(X'*X)^(-1)*(X'*Y);

% applying stepwise regression on Y and X
[a SE PVAL INMODEL STATS NEXTSTEP HISTORY]=stepwisefit(X,Y);
% constant coefficient 
a0_st=STATS.intercept

% plot Y_regr vs Y and provide correlation between them
[r p]=corrcoef([Y Y_regr]);

xline=[min([Y;Y_regr]):max([Y;Y_regr])]; 
figure;
plot(Y,Y_regr,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Y;Y_regr]) max([Y;Y_regr])]);
ylim([min([Y;Y_regr]) max([Y;Y_regr])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['r= ' num2str(r(1,2),'%2.2f')]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2
% Y=a0+a1*X1+a2*X2+a3*X3+a4*X4 + Yrand
%Yrand=randi(100,50,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Yrand.mat;
Ynew=Y+5*Yrand;

% plot Y, Ynew and X
n=size(Y,1);
figure;
plot([1:n],Y,'k-','LineWidth',2);
hold on
plot([1:n],Ynew,'k--','LineWidth',2);
hold on
plot([1:n],X,'LineWidth',1);
legend('Y','Ynew','X1','X2','X3','X4');

% applying MLR on Y and X
b=regress(Ynew,Xnew)
Y_regr=Xnew*b;

% applying stepwise regression on Y and X
[a SE PVAL INMODEL STATS NEXTSTEP HISTORY]=stepwisefit(X,Ynew,'penter',0.05);
% constant coefficient 
a0_st=STATS.intercept

% plot Y_regr vs Y and provide correlation between them
[r p]=corrcoef([Ynew Y_regr]);

xline=[min([Ynew;Y_regr]):max([Ynew;Y_regr])]; 
figure;
subplot(2,1,1)
plot(Ynew,Y_regr,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
ylim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['MLR r= ' num2str(r(1,2),'%2.2f')]);

[ind1 ind2]=find(INMODEL == 1);
b=regress(Ynew,Xnew(:,[1 ind2+1]))
Y_regr=Xnew(:,[1 ind2+1])*b;

% plot Y_regr vs Y and provide correlation between them
[r p]=corrcoef([Ynew Y_regr]);

subplot(2,1,2)
plot(Ynew,Y_regr,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
ylim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['Stepwise regression r= ' num2str(r(1,2),'%2.2f')]);

% check individual correlations (Ynew and each column of X):
[r p]=corrcoef([Ynew X])


% option 2:
% perfroming MLR with all combinations of predictors
% on one part of the data (calibration sample), use another part of the data (validation sample)
% for finding the best model according to R2 (or p-values)
% best model: the one with smallest p-value over the validation sample

% calibration sample: first 25 observations, i.e. Ynew(1:25,:)
% validation sample: the rest, i.e. Ynew(26:end,:)

% set arrey from 1 to the total number of predictors (i.e. 4)
v0 = [1:4];

for kk=1:4
% here pull out all the combinations of 4 numbers (1 to 4) when 1 number is drawn, 
% then 2 numbers (pair) is drawn, then 3 numbers and finally 4 numbers (so this is in the loop)
Ch = nchoosek(v0,kk);

for j=1:length(Ch(:,1))
% calibrate:
B = regress(Ynew(1:25,:),Xnew(1:25,[1 Ch(j,:)+1]));
% validate:
Ytest=Xnew(26:end,[1 Ch(j,:)+1])*B;
[rtest0 ptest0]=corrcoef([Ynew(26:end,:) Ytest]);
rtest(j)=rtest0(1,2);
ptest(j)=ptest0(1,2);
end
[pbest(kk) index(kk)]=min(ptest);
display('combination with columns:')
Ch(index(kk),:)
C(kk).predictors=Ch(index(kk),:);
display('r=') 
Rbest(kk)=rtest(index(kk))
clear ptest
clear rtest
end

% find the best combination and plot the final model
[pbestfinal indexfinal]=min(pbest);  % the best model is the one with the smallest p value
% recalculate the coefficents for this model:
Bfinal = regress(Ynew(1:25,:),Xnew(1:25,[1 C(indexfinal).predictors+1]));
% derive the regressed Y over the whole sample
Yfinal=Xnew(:,[1 C(indexfinal).predictors+1])*Bfinal;

% plot Y_regr vs Y and provide correlation between them
[r p]=corrcoef([Ynew Yfinal]);

figure;
plot(Ynew,Yfinal,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Ynew;Yfinal]) max([Ynew;Yfinal])]);
ylim([min([Ynew;Yfinal]) max([Ynew;Yfinal])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['Optimized regression r= ' num2str(r(1,2),'%2.2f') ' with the predictors ' num2str(C(indexfinal).predictors)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 3
% change coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a0=0;
a1=4;
a2=-3;
a3=2;
a4=-1;

Y=a0+a1*X(:,1)+a2*X(:,2)+a3*X(:,3)+a4*X(:,4);
Ynew_old=Ynew;

Ynew=Y+5*Yrand;

% plot Y, Ynew and X
n=size(Y,1);
figure;
plot([1:n],Ynew_old,'k--','LineWidth',2);
hold on
plot([1:n],Ynew,'r--','LineWidth',2);
hold on
plot([1:n],X,'LineWidth',1);
legend('Ynew old','Ynew','X1','X2','X3','X4');
% applying MLR on Y and X
b=regress(Ynew,Xnew)
Y_regr=Xnew*b;

% applying stepwise regression on Y and X
[a SE PVAL INMODEL STATS NEXTSTEP HISTORY]=stepwisefit(X,Ynew,'penter',0.01);
% constant coefficient 
a0_st=STATS.intercept

% plot Y_regr vs Y and provide correlation between them
[r p]=corrcoef([Ynew Y_regr]);

xline=[min([Ynew;Y_regr]):max([Ynew;Y_regr])]; 
figure;
subplot(2,1,1)
plot(Ynew,Y_regr,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
ylim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['MLR r= ' num2str(r(1,2),'%2.2f')]);

[ind1 ind2]=find(INMODEL == 1);
b=regress(Ynew,Xnew(:,[1 ind2+1]))
Y_regr=Xnew(:,[1 ind2+1])*b;

% plot Y_regr vs Y and provide correlation between them
[r p]=corrcoef([Ynew Y_regr]);

subplot(2,1,2)
plot(Ynew,Y_regr,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
ylim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['Stepwise regression r= ' num2str(r(1,2),'%2.2f')]);

% check individual correlations (Ynew and each column of X):
[r p]=corrcoef([Ynew X])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 4:
% Scale up X4 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X(:,4)=100+10*X(:,4);
Y=a0+a1*X(:,1)+a2*X(:,2)+a3*X(:,3)+a4*X(:,4);

figure;
plot([1:n],Y,'k-','LineWidth',2);
hold on
plot([1:n],X,'LineWidth',1);
legend('Y','X1','X2','X3','X4');

Ynew=Y+5*Yrand;

Xnew=[ones(n,1) X];
% applying MLR on Y and X
b=regress(Ynew,Xnew)
Y_regr=Xnew*b;

% applying stepwise regression on Y and X
[a SE PVAL INMODEL STATS NEXTSTEP HISTORY]=stepwisefit(X,Ynew,'penter',0.01);
% constant coefficient 
a0_st=STATS.intercept

% plot Y_regr vs Y and provide correlation between them
[r p]=corrcoef([Ynew Y_regr]);

xline=[min([Ynew;Y_regr]):max([Ynew;Y_regr])]; 
figure;
subplot(2,1,1)
plot(Ynew,Y_regr,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
ylim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['MLR r= ' num2str(r(1,2),'%2.2f')]);

[ind1 ind2]=find(INMODEL == 1);
b=regress(Ynew,Xnew(:,[1 ind2+1]))
Y_regr=Xnew(:,[1 ind2+1])*b;

% plot Y_regr vs Y and provide correlation between them
[r p]=corrcoef([Ynew Y_regr]);

subplot(2,1,2)
plot(Ynew,Y_regr,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
ylim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['Stepwise regression r= ' num2str(r(1,2),'%2.2f')]);

% standardize X and redo the analysis
% the first column in Xnew consists of ones so the first column should not 
% be rescaled (that is why the loop below goes from 2 to 5)
for k=2:5
Xnew(:,k)=(Xnew(:,k)-mean(Xnew(:,k)))./std(Xnew(:,k));
end

figure;
plot([1:n],Y,'k-','LineWidth',2);
hold on
plot([1:n],Xnew(:,2:end),'LineWidth',1);
legend('Y','X1','X2','X3','X4');

% applying MLR on Y and X
b=regress(Ynew,Xnew)
Y_regr=Xnew*b;

% applying stepwise regression on Y and X
[a SE PVAL INMODEL STATS NEXTSTEP HISTORY]=stepwisefit(Xnew(:,2:end),Ynew,'penter',0.01);
% constant coefficient 
a0_st=STATS.intercept

% plot Y_regr vs Y and provide correlation between them
[r p]=corrcoef([Ynew Y_regr]);

xline=[min([Ynew;Y_regr]):max([Ynew;Y_regr])]; 
figure;
subplot(2,1,1)
plot(Ynew,Y_regr,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
ylim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['MLR r= ' num2str(r(1,2),'%2.2f')]);

[ind1 ind2]=find(INMODEL == 1);
b=regress(Ynew,Xnew(:,[1 ind2+1]))
Y_regr=Xnew(:,[1 ind2+1])*b;

% plot Y_regr vs Y and provide correlation between them
[r p]=corrcoef([Ynew Y_regr]);

subplot(2,1,2)
plot(Ynew,Y_regr,'bo'); hold on
plot(xline,xline,'k-','LineWidth',1);
xlim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
ylim([min([Ynew;Y_regr]) max([Ynew;Y_regr])]);
xlabel('Y');
ylabel('Y_{regr}')
title(['Stepwise regression r= ' num2str(r(1,2),'%2.2f')]);

