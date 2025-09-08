%% Expected value curves under given unc_x (subplot) GAUSSIAN_1D_UNIVARIATE, fig16
% LB,UB,unc_x

LB=-15;
UB=35;
unc_x_list=[0,1.6,3.2,5.8,7,15];
% unc_x_list=[0];

xopt=[-3.00,25.95,20.95,12.00,-10.50]';

for j=1:length(unc_x_list)
    unc_x=unc_x_list(j);
uncx=unc_x/100*(UB-LB)/2;
x_fx=[];
Nsamp=100; %expected value
% Nsamp=2; %edge mean
for xx=LB+uncx:(UB-LB)/1000:UB-uncx
list=linspace(xx-uncx,xx+uncx,Nsamp);
H=[];
for i=1:Nsamp
 H=[H;benchmark_gauss_bimodal_1d(list(i))];
end
x_fx=[x_fx;xx mean(H)];
end

figure(111)
subplot(length(unc_x_list),1,j)
plot(x_fx(:,1),x_fx(:,2))
xlabel('x [-]')
ylabel('𝔼(H(x))')
legend(sprintf('$\\xi_x$=%2.1f\\%%',unc_x),'Interpreter','latex') %reltive in %
legend(sprintf('$|\\xi_x|$=%2.1f',unc_x/100.*(UB-LB)'/2),'Interpreter','latex') %abs.
ylim([0,2])
hold on

if j>1 & j<=length(xopt)
    line(xopt(j-1:j),x_fx(dsearchn(x_fx(:,1),xopt(j-1:j)),2).*ones(2,1),'Color','k','Linestyle','--','Displayname','','Handlevisibility','off','linewidth',0.0005)
end
% x_fx(dsearchn(x_fx(:,1),xopt(j-1:j),2)
end

%saveas(figure(111),'fig16','epsc')
%% Expected value curves under given unc_x (all in one figure) GAUSSIAN_1D_UNIVARIATE
% LB,UB,unc_x

LB=-15;
UB=35;
unc_x_list=[0,1,3,5,10,20];


for j=1:length(unc_x_list)
    unc_x=unc_x_list(j);
uncx=unc_x/100*(UB-LB)/2;
x_fx=[];
Nsamp=100;
for xx=LB+uncx:(UB-LB)/1000:UB-uncx
list=linspace(xx-uncx,xx+uncx,Nsamp);
H=[];
for i=1:Nsamp
 H=[H;benchmark_gauss_bimodal_1d(list(i))];
end
x_fx=[x_fx;xx mean(H)];
end

figure(1)
plot(x_fx(:,1),x_fx(:,2))
xlabel('x [-]')
ylabel('𝔼(H(x))')
lgd{j}=sprintf('$\\xi_x$=%2.0f\\%%',unc_x');
hold on
end
legend(lgd,'Interpreter','latex')


%% Violinplot for evaluating results
folder='Violinplot-Matlab-master'
addpath(genpath(folder))

% load('sensitivity_dx&dw6_violin.mat') %-->results11 %KEINFO

C = colormap('jet');

XXtr=results11{2,7}; %{which  tuning parameter combination,6} %KEINFO
% iterations = 1:1000:length(XXtr); %KEINFO
iterations = [1:10:length(XXtr),length(XXtr)]; %KEINFO

%subplot3
LB=-15;
UB=35;

fx=[];
Nsamp=20; %KEINFO
unc_x=3.2; % percent, %KEINFO (1.6,3.2,1.45abs)
uncx=unc_x/100*(UB-LB)/2; %abs
% uncx=1.15; %abs

forex=XXtr;

figure
for i=1:length(iterations)
    subplot(3,1,1)
vs = Violin({XXtr{iterations(i)}(:,1)},i,'Orientation','vertical','Bandwidth',0.5,'DataStyle','none','Showmedian',false) %'DataStyle', 'histogram','HalfViolin','right');
ylabel('x [-]')
xlim([0.5, length(iterations)+0.5]);
set(gca,xticklabels=num2cell(iterations))
set(gca,xtick=1:length(iterations))
xlabel('Iteration ID [-]');
end
% ylim([min(XXtr{iterations(1)}(:,1)),max(XXtr{iterations(1)}(:,1))])

subplot(3,1,2)
for i=1:length(iterations)
vs = Violin({XXtr{iterations(i)}(:,2)},i,'Orientation','vertical','Bandwidth',0.01,'DataStyle', 'none','Showmedian',false)  %,'DataStyle', 'histogram','HalfViolin','right');
ylabel('H(x) [-]')
xlim([0.5, length(iterations)+0.5]);
set(gca,xtick=1:length(iterations))
set(gca,xticklabels=num2cell(iterations))
xlabel('Iteration ID [-]');
end


xlabel('Iteration ID [-]');
% set correct labels
% set(gca,xticklabels=num2cell(iterations))
% xlim([0.5, length(iterations)+0.5]);

set(gcf,'Units','pixels','Position',[100 100 560 420])





for j=iterations
for i=1:length(forex{j}) %Ns
    xx=forex{j}(i,1); %KEINFO
    list=linspace(xx-uncx,xx+uncx,Nsamp);
    H=[];
    for m=1:Nsamp
     H=[H;benchmark_gauss_bimodal_1d(list(m))];
    end
% fx{j}(i)=benchmark_gauss_bimodal_1d(results11{2,6}{j}(i,1));
fx{j}(i)=mean(H);
end
end


XXtr=fx; %{which  tuning parameter combination,6}

subplot(3,1,3)
for i=1:length(iterations)
% vs = Violin({XXtr{iterations(i)}'},i,'Orientation','vertical','Bandwidth',0.001,'DataStyle', 'histogram','HalfViolin','right');
vs = Violin({XXtr{iterations(i)}'},i,'Orientation','vertical','Bandwidth',0.01,'DataStyle', 'none','Showmedian',false);

ylabel('𝔼(H(x)) [-]')
xlim([0.5, length(iterations)+0.5]);
set(gca,xticklabels=num2cell(iterations))
set(gca,xtick=1:length(iterations))
xlabel('Iteration ID [-]');
end
exportgraphics(figure(2),'fig20.eps','Resolution',300)