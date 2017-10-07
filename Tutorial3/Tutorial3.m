%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 3 (19 Sep 2017)
% PCA on synthetic data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 1
% applying PCA on
% linear progressive wave y(x,t)=sin(k*x-omega*t)
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

% PCA:
[eigenvectors,PCs,eigenvalues]=princomp(y);

% contribution of each mode to total variance:
variance=eigenvalues./sum(eigenvalues);

figure; 
plot([1:100],variance,'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

% plot the first three modes (eigenvestors and PCs)
figure; 
for i=1:3
subplot(3,2,2*i-1)
plot(eigenvectors(:,i));
%ylim([-1 1]);
xlabel('x');
title(['eigenvector',num2str(i)]);

subplot(3,2,2*i)
plot(PCs(:,i));
xlabel('time');
title(['PC',num2str(i)]);
end

% reconstruct data from PCs and eigenvectors
%y(t)-mean(y(t))=sum(PC_j(t)*e_j)

% reconstructing y(1,:) (this is y(t=1))
y_rec(1,:)=PCs(1,1).*eigenvectors(:,1)+PCs(1,2).*eigenvectors(:,2);

% reconstructing all y
for j=2:200
y_rec(j,:)=PCs(j,1).*eigenvectors(:,1)+PCs(j,2).*eigenvectors(:,2);
end

figure;
for i=1:50
subplot(5,10,i)
plot(y_rec(i,:));
ylim([-1 1])
title(['t=',num2str(i)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the same as above but instead using MATLAB function 'princomp'
% apply PCA 'manually' -> as finding eigenvectors and eigenvalues of covariance matrix 
yn=y';

%find mean for each y (1,...m)
mean_y=mean(yn,2);

% subtract each mean from yn
yn=bsxfun(@minus, yn, mean_y);

% calculate covariance
s = cov(yn');

% obtain eigenvalue & eigenvector
[V,D] = eig(s);
eigval = diag(D);
% sort eigenvalues in descending order
eigval = eigval(end:-1:1);

V = fliplr(V);
% get principal components
a=V'*yn;

% contribution of each mode to total variance:
variance2=eigval./sum(eigval);

figure; 
plot([1:100],variance2,'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

% plot the first three modes (eigenvestors and PCs)
figure; 
for i=1:3
subplot(3,2,2*i-1)
plot(V(:,i));
%ylim([-1 1]);
xlabel('x');
title(['e_',num2str(i)]);

subplot(3,2,2*i)
plot(a(i,:));
xlabel('time');
title(['a',num2str(i)]);
end

% compressed way of calculating y_rec from the first two modes
invV=inv(V');
y_rec2=(invV(:,1:2)*a(1:2,:))';

figure;
for i=1:50
subplot(5,10,i)
plot([1:100],y(i,:),'b-',[1:100],y_rec2(i,:),'r-');
if i==1
legend('original','reconstructed');
end
ylim([-1 1])
title(['t=',num2str(i)])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2
% applying PCA on
% linear progressive wave y(x,t)=sin(k*x-omega*t) + noise
% noise=0.5*randn(200,100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load noise.mat

y=y+noise;

figure;
for i=1:50
subplot(5,10,i)
plot(y(i,:));
ylim([-2 2])
title(['t=',num2str(i)])
end

% PCA:
[eigenvectors,PCs,eigenvalues]=princomp(y);

% contribution of each mode to total variance:
variance=eigenvalues./sum(eigenvalues);

figure; 
plot([1:100],variance,'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

% plot the first three modes (eigenvectors and PCs)
figure; 
for i=1:3
subplot(3,2,2*i-1)
plot(eigenvectors(:,i));
xlabel('x');
title(['eigenvector',num2str(i)]);

subplot(3,2,2*i)
plot(PCs(:,i));
xlabel('time');
title(['PC',num2str(i)]);
end

% reconstruct data from PCs and eigenvectors
%y(t)-mean(y(t))=sum(PC_j(t)*e_j)

invE=inv(eigenvectors');
PCf=PCs';
y_rec=(invE(:,1:3)*PCf(1:3,:))';


figure;
for i=1:50
subplot(5,10,i)
plot([1:100],y(i,:),'b-',[1:100],y_rec(i,:),'r-');
ylim([-2 2])
title(['t=',num2str(i)])
if i==1
legend('original','reconstructed');
end
end





