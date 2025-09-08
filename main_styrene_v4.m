%% Validation of the BRANIN objective functions (2024.05.14.)
% close all

% unc=10; % %
% vars={'E_i(2)'}; %ONLY ONE! %KEINFO
% for i=1:length(vars)
%     vars_list=linspace(eval(vars{i})*(1-unc/100),eval(vars{i})*(1+unc/100),4);
% end

abr=1;
vars_list=[1];
LB=[-5,0]';
UB=[5,15]';
LB=[-10,-5]';
UB=[10,20]';


unc_x=[0,0];
points=150; %points to estimate contour plot 30
sum_Z=zeros(points,points);
for j=1:size(vars_list,2)
% eval([vars{1} '=vars_list(j);']);

ment1=[];ment3=[]; 
i=1;
while size(ment1,1)<10000
%Operating point %KEINFO
% xx(1)=-5+ (10-(-5)).*rand(1,1);
% xx(2)=0+ (15-(0)).*rand(1,1);
xx=LB+(UB-LB).*rand(2,1);

% obj=branin(xx + (rand(1,length(LB))*2-1).*unc_x/100.*(LB+UB)'/2); %unc_X
obj=branin(xx);

ment3=[ment3; xx']; %saving the decision variables
ment1=[ment1;obj]; %saving the objective function values
i=i+1;
end

%
%Another plots (03.28.)/a)


if abr==1
%2D: decision variables (P0, SOR)
names1={'x_1','x_2'}; %'SOR [-]'}; %,'T_{steam} [K]'}
[X,Y,Z]=contourZ2(ment3(:,1),ment3(:,2),ment1(:,1),points)
sum_Z=sum_Z+Z;
figg=figure(1)
subplot(2,2,j)
contour(X,Y,Z,points)
xlabel(names1(1))
ylabel(names1(2))
if j==1
legend('Objective function')
end
colorbar;
% title(sprintf('Ei(1)=%6.0f ',E_i_list(j)))
% title([vars{1} sprintf('=%2.0f ',vars_list(j))])
caxis([-30,0])
hold on
end

end

figure
contour(X,Y,sum_Z/size(vars_list,2),points)
xlabel(names1(1))
ylabel(names1(2))
title('Branin function')
colorbar;
caxis([-150,0]);

% figure
% scatter(ment3(:,1),ment3(:,2),[],ment1(:,1))
% colorbar;
% caxis([-100,0]);

%% BRANIN opt. PFO
global unc_x
LB=[-5,0]';
UB=[5,15]';
% unc_x=[3/2.5*100,3/7.5*100,0];
unc_x=[0,0];
tol=[0.001,0.001];
uk_noise=[3 1.25]; %[75/30 1.25]; %[65/30, 1.25] sheel, Ns=50: mean: 121.7164, std: 0.0106; yee: mean:839.95 std:2.12
whub_whlb=300;
Ns=100;
forchange={};
forchange0={};
unc_percent=0;

xopt_robb=[];
for i=1:100
tic
[xopt_rob fxopt M XX XXtr MMM w]=pfo_robust_v2(@branin,Ns,@nonlincon_brannin,LB,UB,'systematic_resampling',0,tol,uk_noise,whub_whlb,forchange,forchange0,unc_percent);
toc
xopt_robb=[xopt_robb;mean(XX(:,1:end-1))]
end
figure;plot(XX(:,1),XX(:,2),'*')
figure;plot(xopt_robb(:,1),xopt_robb(:,2),'*')
xlim([LB(1),UB(1)])
ylim([LB(2),UB(2)])
disp(xopt_robb)
mean(xopt_robb)
std(xopt_robb)
%%
[bb,vv,zz]=uniquetol(xopt_robb,0.5,'Byrows',true);  
[aa,ss,dd]=unique(zz);
figure(10)
h=histogram(zz)
% Get information about the histogram
edges = get(h,'BinEdges');
n = [get(h,'Values'),0]; %above
n1=[0,0]; %below
% Create strings for each bar count
formatSpec = '%.2f';
barstrings = num2str(bb(:,1),formatSpec);
barstrings1 = num2str(bb(:,2),formatSpec);
% Create text objects at each location
% x = (edges(1:end-1)+edges(2:end))/2;
% text(x,n,barstrings,'horizontalalignment','center','verticalalignment','bottom')
% text(x,n,barstrings1,'horizontalalignment','center','verticalalignment','top')
text(x,n1,[['x*=[ ';'x*=[ '],barstrings,[', ';', '],barstrings1,[']';']']],'horizontalalignment','center','verticalalignment','top')
text(x,n,[num2str(n'),['%';'%']],'horizontalalignment','center','verticalalignment','bottom')
set(gca,'xtick',[]);
ylabel('Counts per optimum')
hold on



%% Validation of the GAUSSIAN_BIMODAL objective functions (UNIVARIATE) (2024.07.08.), fig12
close all

abr=1;
vars_list=[1];

%KEINFO
LB=[-15];
UB=[35];


unc_x=[0];
points=20; %points to estimate contour plot
sum_Z=zeros(points,points);

for j=1:size(vars_list,2)
% eval([vars{1} '=vars_list(j);']);

ment1=[];ment3=[]; 
i=1;
while size(ment1,1)<10000
%Operating point %KEINFO
xx=LB+(UB-LB).*rand(1,1);

obj=benchmark_gauss_bimodal_1d(xx + ((rand(1,length(LB))*2-1).*unc_x/100.*(LB+UB)/2)); %unc_X
ment3=[ment3; xx]; %saving the decision variables
ment1=[ment1;obj]; %saving the objective function values
i=i+1;
end
end

figure
plot(ment3,ment1,'.')
xlabel(names1(1))
ylabel('f(x)')
title('Average objective function - expected robust')

%% GAUSSIAN_BIMODAL_1D UNIVARIATE  opt. PFO
global unc_x
LB=-15;
UB=35;

unc_x=[8];
tol=[0.05];
uk_noise=[3 1.25]; %[75/30 1.25]; %[65/30, 1.25] sheel, Ns=50: mean: 121.7164, std: 0.0106; yee: mean:839.95 std:2.12
whub_whlb=10000;
Ns=300;
forchange={};
forchange0={};
unc_percent=0;

xopt_robb=[];
for i=1:20
tic
[xopt_rob fxopt M XX XXtr MMM w]=pfo_robust_v2(@benchmark_gauss_bimodal_1d,Ns,@nonlincon_brannin,LB,UB,'systematic_resampling',1,tol,uk_noise,whub_whlb,forchange,forchange0,unc_percent);

xopt_robb=[xopt_robb;mean(XX(:,1:end-1))]
toc
end
% figure;plot(XX(:,1),XX(:,2),'*')
% uu1=uniquetol(xopt_robb(:,1),0.1);
% histc(xopt_robb,uu1)

disp(xopt_robb)
%
% bb=[0.5,2,9.5,UB]';
% bbb=discretize(xopt_robb,bb)
% h=histogram(bbb)
% % Get information about the histogram
% edges = get(h,'BinEdges');
% n = get(h,'Values');
% % Create strings for each bar count
% barstrings = num2str(bb(1:end-1,1));
% % barstrings1 = num2str(bb(:,2));
% % Create text objects at each location
% x = (edges(1:end-1)+edges(2:end))/2;
% text(x,n,barstrings,'horizontalalignment','center','verticalalignment','bottom')
% % text(x,n,barstrings1,'horizontalalignment','center','verticalalignment','top')
% set(gca,'xtick',[]);
% ylabel('counts per optimum')

[bb,vv,zz]=uniquetol(xopt_robb,0.05,'Byrows',true);  
[aa,ss,dd]=unique(zz);
figure
h=histogram(zz)
% Get information about the histogram
edges = get(h,'BinEdges');
n = get(h,'Values');
% Create strings for each bar count
barstrings = num2str(bb(:,1));
% barstrings1 = num2str(bb(:,2));
% Create text objects at each location
x = (edges(1:end-1)+edges(2:end))/2;
text(x,n,barstrings,'horizontalalignment','center','verticalalignment','bottom')
% text(x,n,barstrings1,'horizontalalignment','center','verticalalignment','top')
set(gca,'xtick',[]);
ylabel('counts per optimum')

figure
subplot(2,1,1)
plot(ment3,ment1,'.')
xlabel('x1')
ylabel('f(x)')
title('Average objective function - expected robust')
subplot(2,1,2)
hist(xopt_robb,50)
xlim([LB, UB])

%% GAUSSIAN_BIMODAL_1D UNIVARIATE  opt. PFO - sensitivity (2024.07.08.); (pfo_robust_v2.m)
%CIKKBE: fig17, fig19, violin...
global unc_x
LB=-15;
UB=35;


tol=[0.05];
uk_noise=[3 1.25]; %[75/30 1.25]; %[65/30, 1.25] sheel, Ns=50: mean: 121.7164, std: 0.0106; yee: mean:839.95 std:2.12
Ns=3000;
forchange={};
forchange0={};
unc_percent=0;


%Sensitivity %KEINFO
% fi_list=[1,1,1,1,0];
% whub_whlb_list=[30000,300,10,3,0];
% unc_x_list=[0,0,0,0,0];
% fi_list=[0,0,0,0,0];
% whub_whlb_list=[0,0,0,0,0];
% unc_x_list=[0,1,3,5,10];

fi_list=[0];
whub_whlb_list=[1500];
unc_x_list=[3.25]; % percent
% fi_list=[1,1,1,1,0,0,0,0,0]; %old strategy; fig19
% whub_whlb_list=[3000,490,100,8,1,1,1,1,1];
% unc_x_list=[0,0,0,0,0,3.25,,4,6];
% fi_list=[1,1,1,1,0]; %old strategy; fig19
% whub_whlb_list=[3000,490,100,8,1];
% unc_x_list=[0,0,0,0,0];

% fi_list=[1,1,1,1,0,0,0,0]; %old strategy
% whub_whlb_list=[30000,300,10,3,0,0,0,0];
% unc_x_list=[0,0,0,0,0,3.5,5,10];
% fi_list=ones(1,5); %new strategy
% whub_whlb_list=3000*ones(1,5);
% unc_x_list=[0,1,3,5,10];

results={};
for sens=1:length(fi_list)
    unc_x=unc_x_list(sens);
xopt_robb=[];
XXtr_save={};
for i=1:20
tic
[xopt_rob fxopt M XX XXtr MMM w]=pfo_robust_v2(@benchmark_gauss_bimodal_1d,Ns,@nonlincon_brannin,LB,UB,'systematic_resampling',fi_list(sens),tol,uk_noise,whub_whlb_list(sens),forchange,forchange0,unc_percent);

xopt_robb=[xopt_robb;mean(XX(:,1:end-1))] %HA uswitsh nincs, MEAN!!
toc

%Concatenate every runs
if isempty(XXtr_save)
    XXtr_save=[XXtr];
else
    JJ=min(length(XXtr_save),length(XXtr));
    for j=1:JJ %iterations
    XXtr_JJ=XXtr_save{JJ};
    XXtr_save{j}=[XXtr_save{j};XXtr{j}];
    end
        for j=JJ+1:max(length(XXtr_save),length(XXtr))
            if length(XXtr_save)<length(XXtr)
%             XXtr_save{j}=[XXtr{j}];
            XXtr_save{j}=[XXtr_JJ;XXtr{j}];
            elseif length(XXtr_save)>length(XXtr)
            XXtr_save{j}=[XXtr_save{j};XXtr{JJ}];
            end
    end
end
end
% figure;plot(XX(:,1),XX(:,2),'*')
% uu1=uniquetol(xopt_robb(:,1),0.1);
% histc(xopt_robb,uu1)

disp(xopt_robb)
results{sens,1}=xopt_robb;
results{sens,2}=fi_list(sens);
results{sens,3}=whub_whlb_list(sens);
results{sens,4}=unc_x;

[bb,vv,zz]=uniquetol(xopt_robb,0.05,'Byrows',true);  
[aa,ss,dd]=unique(zz);
figure
h=histogram(zz)
% Get information about the histogram
edges = get(h,'BinEdges');
n = get(h,'Values');
% Create strings for each bar count
barstrings = num2str(bb(:,1));
% barstrings1 = num2str(bb(:,2));
% Create text objects at each location
x = (edges(1:end-1)+edges(2:end))/2;
text(x,n,barstrings,'horizontalalignment','center','verticalalignment','bottom')
% text(x,n,barstrings1,'horizontalalignment','center','verticalalignment','top')
set(gca,'xtick',[]);
ylabel('counts per optimum')

results{sens,5}=Ns;
results{sens,6}=XXtr; %all particles (x,fx) in all iteration
results{sens,7}=XXtr_save; %all particles (x,fx) in all iteration in all 20 runs

% disp(XXtr_save{end})
% figure;hist(XXtr{end}(:,1));hist(XXtr_save{end}(:,1))
end



%% FIGURE fig19 (sensitivity_dx&dw1.mat)
load('sensitivity_dx&dw1.mat')
results=results11;

yyy=[];
xxx=[];
edges=[-20,-5,3,15,22.5,30];

% load('sensitivity_Ux.mat')
% for i=1:length(results)
% yyy=[yyy;histc(results{i,1}, edges)'];
% xxx=[xxx;results{i,4}];
% end

% load('sensitivity_Ux&dw.mat')
for i=1:length(results)
yyy=[yyy;histc(results{i,1}, edges)'];
xxx=[xxx;results{i,4}-results{i,3}];
end



figure(19)
% bar(xxx,yyy)
bar(yyy)
xlabel('($\Delta w$, $\Delta x$)','Interpreter','latex')
ylabel('Counts [-]')
% title('N_s=3000')
% set(gca, 'XTickLabel', {'Model1' 'Model2'})
% xticklabels({'(3000,0)','(490,0)','(100,0)','(8,0)','(1,0)','(1,3.25)','(1,3.7)','(1,4)','(1,6)'})
% %rel. kszi
% xticklabels({'(3000, 0)','(490, 0)','(100, 0)','(8, 0)','(1, 0)','(1, 0.81)','(1, 0.93)','(1, 1)','(1, 1.5)'}) %abs kszi
xticklabels({'(3000, 0)','(490, 0)','(100, 0)','(8, 0)','(1, 0)','(1, 0.8)','(1, 0.9)','(1, 1.0)','(1, 1.5)'}) %abs kszi


% legend('-10.55','-3.02','12.06','20.95','25.94') %opt
legg=legend('-10.50','-3.00','12.00','20.95','25.95') %cfgv
title(legg,'x^*_{robust}=','interpreter','tex')

% saveas(figure(19),'fig19','epsc')
%%
[~,ind]=sort(ment3);
results=results([results{:,4}]~=3,:);

figure
subplot(size(results,1)+3,1,1:2)
% plot(ment3,ment1,'.')
plot(ment3(ind),ment1(ind))
xlabel({'x'})
ylabel('H(x)')
hold on
for i=1:size(results,1)
subplot(size(results,1)+3,1,i+3)
hist(results{i,1},linspace(LB,UB,50))
xlim([LB, UB])
ylim([0,100])
xlabel({'x*'})
ylabel('Counts') % per optimum')
if results{i,2}==1
    if i==1
    title({sprintf('\\sf Exponential  $\\varphi_k$ ($\\Delta w$=%1.0f)',results{i,3})},'Interpreter','latex');
%     title(sprintf('Exponential $\varphi_k$ (%2.0f)',results{i,2}),'Interpreter','latex')
    end
        title({sprintf('\\sf Exponential  $\\varphi_k$ ($\\Delta w$=%1.0f)',results{i,3})},'Interpreter','latex');
elseif results{i,2}==0
    title({sprintf('\\sf Uniform  $\\varphi_k$ ($U_x$=%1.0f\\%%)',results{i,4})},'Interpreter','latex');
%     title('\sf \textbf{Uniform} $\varphi_k$','Interpreter','latex')
end
end



%% GAUSSIAN_BIMODAL_1D UNIVARIATE  opt. PFO - Ns sensitivity??? (2024.07.15.)
global unc_x
LB=-15;
UB=35;


tol=[0.05];
uk_noise=[3 1.25]; %[75/30 1.25]; %[65/30, 1.25] sheel, Ns=50: mean: 121.7164, std: 0.0106; yee: mean:839.95 std:2.12
Ns_list=[300];
forchange={};
forchange0={};
unc_percent=0;


%Sensitivity %KEINFO
% fi_list=[1,1,1,1,0];
% whub_whlb_list=[30000,300,10,3,0];
% unc_x_list=[0,0,0,0,0];
% fi_list=[0,0,0,0,0];
% whub_whlb_list=[0,0,0,0,0];
% unc_x_list=[0,1,3,5,10];
fi_list=[1];
whub_whlb_list=[10];
unc_x_list=[0];

% fi_list=[1,1,1,1,0,0,0,0]; %old strategy CIKK
% whub_whlb_list=[30000,300,10,3,0,0,0,0];
% unc_x_list=[0,0,0,0,0,3,5,10];

% fi_list=ones(1,5); %new strategy
% whub_whlb_list=3000*ones(1,5);
% unc_x_list=[0,1,3,5,10];

results={};
for sens1=1:length(fi_list)
    unc_x=unc_x_list(sens1);
    ddx=unc_x/100.*(UB-LB)'/2;
for sens=1:length(Ns_list)
    Ns=Ns_list(sens);
xopt_robb=[];
for i=1:20
tic
[xopt_rob fxopt M XX XXtr MMM w]=pfo_robust_v2(@benchmark_gauss_bimodal_1d,Ns,@nonlincon_brannin,LB,UB,'systematic_resampling',fi_list(sens1),tol,uk_noise,whub_whlb_list(sens1),forchange,forchange0,unc_percent);

xopt_robb=[xopt_robb;mean(XX(:,1:end-1))] %HA uswitsh nincs, MEAN!!
toc
end
% figure;plot(XX(:,1),XX(:,2),'*')
% uu1=uniquetol(xopt_robb(:,1),0.1);
% histc(xopt_robb,uu1)

disp(xopt_robb)
results{sens,1}=xopt_robb;
results{sens,2}=fi_list(sens1);
results{sens,3}=whub_whlb_list(sens1);
results{sens,4}=unc_x;

[bb,vv,zz]=uniquetol(xopt_robb,0.05,'Byrows',true);  
[aa,ss,dd]=unique(zz);
figure
h=histogram(zz)
% Get information about the histogram
edges = get(h,'BinEdges');
n = get(h,'Values');
% Create strings for each bar count
barstrings = num2str(bb(:,1));
% barstrings1 = num2str(bb(:,2));
% Create text objects at each location
x = (edges(1:end-1)+edges(2:end))/2;
text(x,n,barstrings,'horizontalalignment','center','verticalalignment','bottom')
% text(x,n,barstrings1,'horizontalalignment','center','verticalalignment','top')
set(gca,'xtick',[]);
ylabel('counts per optimum')

results{sens,5}=Ns;
results{sens,6}=max(n);
end
end

%%
[~,ind]=sort(ment3);
results=results([results{:,4}]~=3,:);

figure
subplot(size(results,1)+3,1,1:2)
% plot(ment3,ment1,'.')
plot(ment3(ind),ment1(ind))
xlabel({'x'})
ylabel('H(x)')
hold on
for i=1:size(results,1)
subplot(size(results,1)+3,1,i+3)
hist(results{i,1},linspace(LB,UB,50))
xlim([LB, UB])
ylim([0,100])
xlabel({'x*'})
ylabel('Counts') % per optimum')
if results{i,2}==1
    if i==1
    title({sprintf('\\sf Exponential  $\\varphi_k$ ($\\Delta w$=%1.0f)',results{i,3})},'Interpreter','latex');
%     title(sprintf('Exponential $\varphi_k$ (%2.0f)',results{i,2}),'Interpreter','latex')
    end
        title({sprintf('\\sf Exponential  $\\varphi_k$ ($\\Delta w$=%1.0f)',results{i,3})},'Interpreter','latex');
elseif results{i,2}==0
    title({sprintf('\\sf Uniform  $\\varphi_k$ ($U_x$=%1.0f\\%%)',results{i,4})},'Interpreter','latex');
%     title('\sf \textbf{Uniform} $\varphi_k$','Interpreter','latex')
end
end






%% Functions

function [c,ceq] = nonlincon_brannin(x) %=boundary(x)
        c=-1*zeros(1,2);
        ceq = [];
        %c=[];
        
end

