figure;ksdensity(fx{k}(X{k}(:,1)<1.97));hold on;ksdensity(fx{k}(X{k}(:,1)>1.97))
title('k=47'),xlabel('fx');legend('x1<1.97','x1>1.97')

%Objective value of the particles by x1, x2...
figure;subplot(2,1,1);scatter(X{k}(:,1),fx{k});xlabel('x1'),ylabel('fx'),title(sprintf('k=%2.0f',k))
subplot(2,1,2);scatter(X{k}(:,2),fx{k});xlabel('x2'),ylabel('fx')

%Objective value of the particles by x1, x2... 
figure;subplot(3,1,1);scatter(X{k}(:,1),fx{k});xlabel('x1'),ylabel('fx'),title(sprintf('k=%2.0f',k))
subplot(3,1,2);scatter(X{k}(:,2),fx{k});xlabel('x2'),ylabel('fx');
subplot(3,1,3);scatter(X{k}(:,3),fx{k});xlabel('x3'),ylabel('fx');

%3D plot of objective values by (x1,x2)
X=XXtr;
fx={};
fx{k}=XXtr{k}(:,end);
figure;scatter3(X{k}(:,1),X{k}(:,2),fx{k},[],fx{k},'filled');xlabel('x1'),ylabel('x2'),title(sprintf('k=%2.0f',k));col=colorbar;ylabel(col,'fx');caxis([50*10.55,200*10.55])


%branin
X=XXtr;
fx={};
fx{k}=XXtr{k}(:,end);
figure;scatter3(X{k}(:,1),X{k}(:,2),fx{k},[],fx{k},'filled');xlabel('x1'),ylabel('x2'),title(sprintf('k=%2.0f',k));col=colorbar;ylabel(col,'fx');caxis([-5,0])


%Measure of exploration (u_switch threshold)
(max(X{k}(:,1))-min(X{k}(:,1)))/tol(1)*(max(X{k}(:,2))-min(X{k}(:,2)))/tol(2)

%Post: (exploitation): number of particles in the certain points
[uu1]=unique(X{k}(:,1)); histc(X{k},uu1)

[uu1,uu2,uu3]=unique(X{k},'rows');
numel(find(uu3==3))


% uniquetol(X{k},0.5,','ByRows',true)

%%
clust = zeros(size(X{k},1),6);
for i=1:6; tic;
clust(:,i) = kmeans(X{k},i)%,'emptyaction','singleton','replicate',5);
end;toc
va = evalclusters(X{k},clust,'CalinskiHarabasz')

tic;eva=evalclusters(X{k},'kmeans','gap','KList',1:6);toc