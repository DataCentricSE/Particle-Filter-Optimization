
options = optimoptions('particleswarm','SwarmSize',100)

Xopt=[];
for i=1:100
X=particleswarm(@branin,2,[-5,0],[5,15],options);
Xopt=[Xopt;X]
% X=fmincon(@branin,[3,7.5],[],[],[],[],[-5,0],[5,15])
end

uu1=uniquetol(Xopt(:,1),0.1);
hist(Xopt(:,1),uu1)

%%
[bb,vv,zz]=uniquetol(Xopt,0.5,'Byrows',true);  
[aa,ss,dd]=unique(zz);
figure(10)
h=histogram(zz)
% Get information about the histogram
edges = get(h,'BinEdges');
n = get(h,'Values'); %above
n1=[0,0]; %below
% Create strings for each bar count
formatSpec = '%.2f';
barstrings = num2str(bb(:,1),formatSpec);
barstrings1 = num2str(bb(:,2),formatSpec);
% Create text objects at each location
x = (edges(1:end-1)+edges(2:end))/2;
% text(x,n,barstrings,'horizontalalignment','center','verticalalignment','bottom')
% text(x,n,barstrings1,'horizontalalignment','center','verticalalignment','top')
text(x,n1,[['x*=[ ';'x*=[ '],barstrings,[', ';', '],barstrings1,[']';']']],'horizontalalignment','center','verticalalignment','top')
text(x,n,[num2str(n'),['%';'%']],'horizontalalignment','center','verticalalignment','bottom')
set(gca,'xtick',[]);
ylabel('Counts per optimum')
hold on
legend('PFO','PSO')
bins=unique(zz);

%Complete with PFO plot
%%
[bb,vv,zz]=uniquetol(xopt_robb,0.5,'Byrows',true);  
[aa,ss,dd]=unique(zz);
figure(10)
h=histogram(zz*2)
% Get information about the histogram
edges = get(h,'BinEdges');
n = [0,get(h,'Values')]; %above
n1=[0,0]; %below
% Create strings for each bar count
formatSpec = '%.2f';
barstrings = num2str([0,mean(xopt_robb(zz==1,1))],formatSpec);
barstrings1 = num2str([0,mean(xopt_robb(zz==1,2))],formatSpec);
% Create text objects at each location
% x = (edges(1:end-1)+edges(2:end))/2;
% text(x,n,barstrings,'horizontalalignment','center','verticalalignment','bottom')
% text(x,n,barstrings1,'horizontalalignment','center','verticalalignment','top')
% text(x,n1,[['x*=[ ';'x*=[ '],barstrings,[', ';', '],barstrings1,[']';']']],'horizontalalignment','center','verticalalignment','top')
text(x,n,[num2str(n'),['%';'%']],'horizontalalignment','center','verticalalignment','bottom')
set(gca,'xtick',[]);
ylabel('Counts per optimum')
hold on

function yy=branin(x)
a=1;
b=5.1/(4*pi()^2);
c=5/pi();
r=6;
s=10;
t=1/8/pi();
yy=(a*(x(2)-b*x(1)^2+c*x(1)-r)^2+s*(1-t)*cos(x(1))+s);
end