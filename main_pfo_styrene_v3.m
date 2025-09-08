clear all
close all
clc

global G0 d Nc Nr L rho_b D_p e A_t nu_ij E_i A_i M_j G a_i b_i alfa_j beta_j gamma_j a_j b_j c_j ment ment1
global F_st0 F_be0 F_to0 abr
global H_st H_eb H_be H_tol H_steam
global C_st C_be C_to C_h1 C_h2 ment2 ment3 abr ment7
global G_st G_to G_be G_H2 G_meth G_eth G_steam1 G_steam2 E_i A_i
global LB UB unc_x

% Available data
%j=[ethylbenzene (eb),styrene (st),benzene (be),toluene (tol),ethylene (eth),methane (meth),steam,H2,CO,CO2] 
%i=1,...,Nc

abr=1;
%Parameters of the reaction system
Nc=10; %number of components
Nr=6; %number of reactions
nu_ij=[-1 1 0 0 0 0 0 1 0 0;
       -1 0 1 0 1 0 0 0 0 0;
       -1 0 0 1 0 1 0 -1 0 0;
       [0 0 0 0 -1 0 -2 4 2 0]/2; %elnashaie1994 p.195
       0 0 0 0 0 -1 -1 3 1 0;
       0 0 0 0 0 0 -1 1 -1 1];
M_j=[106.167 104.152 78.114 92.141 28.054 16.043 18.015 2.016 28.01 44.01]; %kg/kmol, molar mass

%Input reactor specification
d=1.95; %m, reactor diameter
L=1.7; %m, reactor length
rho_b=2146; %kg/m3, catalyst bulk density
D_p=0.0047; %m, catalyst particle diameter
e=0.445; %-, bed-void fraction

%Inlet specification
T0=922.59; %K
P0=2.4; %bar
F_eb0=36.87/3600; %kmol/s, ethyl-benzene

F_steam0=453.1/3600; %kmol/s, steam

%Calculated data
A_t=d^2*pi()/4;

% Specific heat (organic: cp_j=alfa_j+beta_j*T+gamma_j*T^2; inorganic: cp_j=a_j+b_j*T+c_j/T^2)
%j=[ethylbenzene (eb),styrene (st),benzene (be),toluene (tol),ethylene (eth),methane (meth),steam,H2,CO,CO2] 

%smihtvanness1975
alfa_j=[9.34 17.04 -1.71 2.41 11.84 14.15]; %kJ/kmol/K
beta_j=[0.460 0.417 0.325 0.391 0.120 0.075]; %kJ/kmol/K^2
gamma_j=-10^(-5)*[15.36 13.85 11.06 13.07 3.65 1.8]; %kJ/kmol/K^3
% 
a_j=[28.85 27.01 28.07 45.37]; %kJ/kmol/K
b_j=[0.01206 3.51e-3 4.63e-3 8.69e-3]; %kJ/kmol/K^2
c_j=[1.01 0.69 -0.26 -9.62]*10^5; %kJ*K/kmol

%Reaction heat
a_i=[120750.11 108818.11 -53178.2 82065.74 211255.19 -45223.66]; %elsőnél átváltás? (sheel1969,elnashaie1994)
b_i=[4.56 -7.96 -13.19 8.84 16.58 10.47];

%Reaction rate constants
%yee2003,chaudhari2022
E_i=[90981.4 207989.2 91515.3 103996.7 65723.3 73628.4]; %kJ/kmol
A_i=[-0.0854 13.2392 0.2961 -0.0724 -2.9344 21.2402];

%yee2003 profit
%Parameters
H_st=103; %$/kmol
H_be=30.8; %$/kmol
H_tol=33.9; %$/kmol
H_eb=45.6; %$/kmol
H_steam=0.36; %$/kmol

%sheel1969 profit
%Parameters ($/kmol)
G_st=10.10;
G_to=-10.36;
G_be=-7.75;
G_H2=0.13;
G_meth=0.42;
G_eth=0.68;
G_steam1=-4.69e-3;
G_steam2=-2.57e-5;

%% fmincon (adiabatic)
if abr==0
LB=[1 7 850];
UB=[2.63 20 1200]; 
x0_fm=((LB+UB)/2-LB)./(UB-LB);
tic
[x_opt_fm fval_fm exitflag_fm]=fmincon(@optsystem_foropt_fm,x0_fm,[],[],[],[],zeros(3,1)',ones(3,1)',@nonlincon_fm)
toc
x_opt_fm=x_opt_fm.*(UB-LB)+LB
end


%% Sensitivity analysis of the optimum (bifurcation analysis, fmincon, adiabatic,obj3)
global LB UB E_i A_i
global C_st C_be C_to C_h1 C_h2
global G_st G_be G_to G_steam2 G_H2 G_meth G_eth G_steam1

if abr==0
Nx=3;
Np=7;
ment4=cell(1,Np);ment8=cell(1,Np);
Npoints=10;
Nruns=1;

LB=[1 7 850];
UB=[2.63 20 1200];  %KE: INNEN
x0_fm=((LB+UB)/2-LB)./(UB-LB); %normalized

for j=1:Np
for i=1:Npoints
    G_st=10.10*10.5518; %INFLATED
    G_to=-10.36*10.5518;
    G_be=-7.75*10.5518;
    G_H2=0.13*10.5518;
    G_meth=0.42*10.5518;
    G_eth=0.68*10.5518;
    G_steam1=-4.69e-3*10.5518;
    G_steam2=-2.57e-5*10.5518;


        G_st=123.84;
    G_to=-85.04;
    G_be=-104.39;
    G_H2=0.13;
    G_meth=0.42;
    G_eth=0.68;
    G_steam1=-10*4.69e-3;
    G_steam2=-10*2.57e-5;


    E_i=[90981.4 207989.2 91515.3 103996.7 65723.3 73628.4]; %kJ/kmol
    A_i=[-0.0854 13.2392 0.2961 -0.0724 -2.9344 21.2402];
   if ismember(j,1:4)
        elt=0.15;
    else 
        elt=0.10;
    end
 
    HHH=[G_st G_be G_to G_steam2 E_i(1) E_i(2) A_i(2)];
    HHH_list=linspace(HHH(j)*(1-elt),HHH(j)*(1+elt),Npoints);
    HHH(j)=HHH_list(i);
    HH=num2cell(HHH);
    
    %plus=(rand(1,5)-0.5*ones(1,5)).*(2*[H_st H_be H_tol H_eb H_steam]*elt);
    %plus(1:end~=j)=0; %changing only one parameter
    %HH=num2cell([H_st H_be H_tol H_eb H_steam]+plus);
    
    [G_st G_be G_to G_steam2 E_i(1) E_i(2) A_i(2)]=deal(HH{:}); %randomize the parameters in [par-par/elt par+par/elt]
    ment4{j}=[ment4{j}; G_st G_be G_to G_steam2*1000 E_i(1)/1000 E_i(2)/1000 A_i(2)]; %saving the profit parameters (MJ és 10^-3 $ Gsteam2-nél!!!)
    %27.56
    
    ment9=[];
    for ii=1:Nruns
        % x0_fm=rand(Nx,1);
        tic
        [x_opt_fm fval_fm exitflag_fm]=fmincon(@optsystem_foropt_fm,x0_fm,[],[],[],[],zeros(Nx,1)',ones(Nx,1)',@nonlincon_fm);
        toc;
        % x0_fm=x_opt_fm;
        x_opt_fm=(x_opt_fm.*(UB-LB)+LB).*[1 1 1];
        ment9=[ment9; x_opt_fm fval_fm exitflag_fm toc];
    end

    [~,idx]=min(ment9(:,end-1));  
    ment8{j}=[ment8{j}; ment9(idx,:)]
end
end
end

%% 
abr=0
if abr==0
names1={'P_0 [bar]','SOR [-]','T_{steam} [K]', 'H(x) [USD/h]'};
% names3={'H_s_t [$/kmol]', 'H_b_e [$/kmol]','H_t_o_l [$/kmol]','H_e_b [$/kmol]','H_{steam} [$/kmol]','H_{energy} [$/kJ]'};
%names3={'G_s_t [$/kmol]', 'G_b_e [$/kmol]','G_t_o_l [$/kmol]','G_{steam1} [$/kmol]','G_{steam2} [$/kmol]','H_{energy} [$/kJ]'};
names3={'G_s_t [$/kmol]', 'G_b_e [$/kmol]','G_t_o [$/kmol]','G_{steam2} [10^{-3}$/kmol]','E_i(1) [MJ/kmol]','E_i(2) [MJ/kmol]','A_i(2)'  }

LB=[1 7 850 0];
UB=[2.63 20 1200 2000];
limits=[LB' UB'];


figure
for i=1:Np %# of the changed parameter
for j=1:Nx+1 %the decision variables in the optimum
    subplot(size(ment4,2),size(x_opt_fm,2),(i-1)*size(x_opt_fm,2)+j)
    plot(ment4{i}(:,i),ment8{i}(:,j),'.')
    xlabel(names3{i})
    ylabel(names1{j})
    ylim(limits(j,:))
    hold on
end
end
end

%% FIGURE 1 : 
clear all; load('save_sheel_bifurcation_PFO.mat')
abr=0
if abr==0
names1={'P_0 [bar]','SOR [-]','T_{steam} [K]', 'H(x) [USD/h]'};
% names3={'H_s_t [$/kmol]', 'H_b_e [$/kmol]','H_t_o_l [$/kmol]','H_e_b [$/kmol]','H_{steam} [$/kmol]','H_{energy} [$/kJ]'};
%names3={'G_s_t [$/kmol]', 'G_b_e [$/kmol]','G_t_o_l [$/kmol]','G_{steam1} [$/kmol]','G_{steam2} [$/kmol]','H_{energy} [$/kJ]'};
names3={'G_s_t [$/kmol]', 'G_b_e [$/kmol]','G_t_o [$/kmol]','G_{steam2} [10^{-3}$/kmol]','E_i(1) [MJ/kmol]','E_i(2) [MJ/kmol]','A_i(2) [-]'  }

LB=[1 7 850 0];
UB=[2.63 20 1200 2000];
limits=[LB' UB'];

x_opt_fm=zeros(1,5);

figure
for i=1:Np %# of the changed parameter
    ment8{i}(:,end)=ment8{i}(:,end)*10.5518; %inflation
    ment4{i}(:,1:4)=ment4{i}(:,1:4)*10.5518; %inflation
    ment4{i}(:,4)=ment4{i}(:,4)*1000; %1e-3
    ment4{i}(:,5:6)=ment4{i}(:,5:6)/1000; %kJ
for j=1:Nx+1 %the decision variables in the optimum
    subplot(size(ment4,2),size(x_opt_fm,2)-1,(i-1)*(size(x_opt_fm,2)-1)+j)
    plot(ment4{i}(:,i),ment8{i}(:,j),'.')
    xlabel(names3{i})
    ylabel(names1{j})
    ylim(limits(j,:))
    hold on
end
end
end
%% PFO test on the styrene reactor
%x is not normalized (az algoritmus jellegéből fakadóan itt mindegy)
if abr==0
%KE: hangolandók:
tol=0.05;
uk_noise=[75/30 1.25]; %[65/30, 1.25] sheel, Ns=50: mean: 121.7164, std: 0.0106; yee: mean:839.95 std:2.12
whub_whlb=300;
Ns=150; 


LB=[1 7 850]'; 
UB=[2.63 20 1200]'; 
perf=[];
Xopt=[];
abr=1;

%bbb_list=40:5:90;

% for j=1:length(bbb_list)
%     bbb=bbb_list(j);
%     perf=[];
for i=1:10
    tic
[xopt fxopt M XX]=pfo(@optsystem_foropt_pfo,Ns,@nonlincon_pfo,LB,UB,'systematic_resampling',1,tol,uk_noise,whub_whlb);
toc 
perf=[perf fxopt]
Xopt=[Xopt; xopt];
end

mean(perf)
std(perf)
%PP(j)=mean(perf);
%SS(j)=std(perf);
% end

% figure;plot(XX{4}(:,5),'.');ylabel('x_5');xlabel('# of the particle')
end
%% PFO bifurcation analysis
if abr==1
global LB UB
global H_st H_eb H_be H_tol H_steam H_energy

%KE: hangolandók:
tol=0.05;
uk_noise=[75/30 1.25];% 140/30 2.3];
whub_whlb=300;
Ns=150; 
Np=7;

LB=[1 7 850]';
UB=[2.63 20 1200]'; %KE: INNEN

perf=[];
Xopt=[];
abr=1;

ment4=cell(1,Np);ment8=cell(1,Np);
Npoints=10;

%Parameters
for j=1:Np
for i=1:Npoints
    G_st=10.10*10.5518;
    G_to=-10.36*10.5518;
    G_be=-7.75*10.5518;
    G_H2=0.13*10.5518;
    G_meth=0.42*10.5518;
    G_eth=0.68*10.5518;
    G_steam1=-4.69e-3*10.5518;
    G_steam2=-2.57e-5*10.5518;
    E_i=[90981.4 207989.2 91515.3 103996.7 65723.3 73628.4]; %kJ/kmol
    A_i=[-0.0854 13.2392 0.2961 -0.0724 -2.9344 21.2402];
% 
% H_st=103; %$/kmol
% H_be=30.8; %$/kmol
% H_tol=33.9; %$/kmol
% H_eb=45.6; %$/kmol
% H_steam=0.36; %$/kmol

if ismember(j,1:4)
    elt=0.15;
else 
    elt=0.10;
end

% HHH=[H_st H_be H_tol H_eb H_steam];
HHH=[G_st G_be G_to G_steam2 E_i(1) E_i(2) A_i(2)]; %[G_st G_be G_to G_steam1 G_steam2 A_i(2) E_i(1) G_meth G_eth G_H2];
HHH_list=linspace(HHH(j)*(1-elt),HHH(j)*(1+elt),Npoints);
HHH(j)=HHH_list(i);
HH=num2cell(HHH);

[G_st G_be G_to G_steam2 E_i(1) E_i(2) A_i(2)]=deal(HH{:}); %randomize the parameters in [par-par/elt par+par/elt]
ment4{j}=[ment4{j}; G_st G_be G_to G_steam2 E_i(1) E_i(2) A_i(2)]; %saving the profit parameters
%27.56

tic
[xopt fxopt M XX]=pfo(@optsystem_foropt_pfo,Ns,@nonlincon_pfo,LB,UB,'systematic_resampling',1,tol,uk_noise,whub_whlb);
tt(j,i)=toc
ment8{j}=[ment8{j}; xopt fxopt]
end
end
end

%%
if abr==0
names1={'P_0 [bar]','SOR [-]','T_{steam} [K]','delta [%]','lambda [%]'};
% names3={'H_s_t [$/kmol]', 'H_b_e [$/kmol]','H_t_o_l [$/kmol]','H_e_b [$/kmol]','H_{steam} [$/kmol]','H_{energy} [$/kJ]'};
%names3={'G_s_t [$/kmol]', 'G_b_e [$/kmol]','G_t_o_l [$/kmol]','G_{steam1} [$/kmol]','G_{steam2} [$/kmol]','A_i(1)', 'E_i(1) [kJ/kmol]' ,'G_{meth} [$/kmol]' ,'G_{eth} [$/kmol]', 'G_{H2} [$/kmol]'};
names3={'G_s_t [$/kmol]', 'G_b_e [$/kmol]','G_t_o [$/kmol]','G_{steam2} [10^{-3}$/kmol]','E_i(1) [MJ/kmol]','E_i(2) [MJ/kmol]','A_i(2)'  }


LB=[1 7 850 0]';
UB=[2.63 20 1200 200]';
limits=[LB UB];


figure
for i=1:Np %# of the changed parameter
for j=1:Nx %the decision variables in the optimum
    subplot(Np,Nx,(i-1)*Nx+j)
    plot(ment4{i}(:,i),ment8{i}(:,j),'.')
    xlabel(names3{i})
    ylabel(names1{j})
    ylim(limits(j,:))
    hold on
end
end

if abr==0
figure
for i=size(ment4,2)/2+1:size(ment4,2) %# of the changed parameter
for j=1:size(xopt,2) %the decision variables in the optimum
    subplot(size(ment4,2)/2,size(xopt,2),(i-6)*size(xopt,2)+j)
    plot(ment4{i}(:,i),ment8{i}(:,j),'.')
    xlabel(names3{i})
    ylabel(names1{j})
    ylim(limits(j,:))
    hold on
end
end
end
end

figure
for i=1:Np %# of the changed parameter
for j=1:Nx+1 %the decision variables in the optimum
    subplot(size(ment4,2),size(x_opt_fm,2),(i-1)*size(x_opt_fm,2)+j)
    plot(ment4{i}(:,i),ment8{i}(:,j),'.')
    xlabel(names3{i})
    ylabel(names1{j})
    ylim(limits(j,:))
    hold on
end
end

%% Robust PFO (pfo_robust_v1.m)
%Set back data
G_st=10.10;
G_to=-10.36;
G_be=-7.75;
G_H2=0.13;
G_meth=0.42;
G_eth=0.68;
G_steam1=-4.69e-3;
G_steam2=-2.57e-5;

E_i=[90981.4 207989.2 91515.3 103996.7 65723.3 73628.4]; %kJ/kmol
A_i=[-0.0854 13.2392 0.2961 -0.0724 -2.9344 21.2402];



if abr==1
%KE: hangolandók:
tol=[0.02,0.1,2];
uk_noise=[3 1.25]; %[75/30 1.25]; %[65/30, 1.25] sheel, Ns=50: mean: 121.7164, std: 0.0106; yee: mean:839.95 std:2.12
whub_whlb=300;
Ns=250; 


forchange={'E_i(2)','E_i(1)','A_i(2)','G_st'};  %KE:INFO
unc_percent=[0,10,0,10]; % in % %KE:INFO


forchange0=[];
for i=1:length(forchange)
    forchange0=[forchange0,eval(forchange{i})]; %original/predicted/supposed value
end



LB=[1 7 1037.3]'; 
UB=[2.63 20 1037.3]'; 
perf=[];
Xopt=[];

for i=1:1
    tic
[xopt_rob fxopt M XX XXtr MMM w]=pfo_robust_v1(@optsystem_foropt_pfo,Ns,@nonlincon_pfo,LB,UB,'systematic_resampling',0,tol,uk_noise,whub_whlb,forchange,forchange0,unc_percent);
toc 
perf=[perf fxopt]
Xopt=[Xopt; xopt_rob];
end

% mean(perf)
% std(perf)
end
xopt_robb=[mean(XX(:,1)) mean(XX(:,2)) mean(XX(:,3))]

%% Robust PFO (pfo_robust_v2.m; with uncertainty on x, 2024.05.09.)
%Set back data
G_st=10.10*10.5518; %$/kmol
G_to=-10.36*10.5518;
G_be=-7.75*10.5518;
G_H2=0.13*10.5518;
G_meth=0.42*10.5518;
G_eth=0.68*10.5518;
G_steam1=-4.69e-3*10.5518;
G_steam2=-2.57e-5*10.5518;

E_i=[90981.4 207989.2 91515.3 103996.7 65723.3 73628.4]; %kJ/kmol
A_i=[-0.0854 13.2392 0.2961 -0.0724 -2.9344 21.2402];



if abr==1
%KE: hangolandók:
tol=[0.05,0.1,2];
uk_noise=[3 1.25]; %[75/30 1.25]; %[65/30, 1.25] sheel, Ns=50: mean: 121.7164, std: 0.0106; yee: mean:839.95 std:2.12
whub_whlb=300;
Ns=150; 


forchange={'E_i(2)','E_i(1)','A_i(2)','G_st','G_to','G_be','G_steam2','E_i(3)'};  %KE:INFO
unc_percent=[0,0,0,0,0,0,0,0]; % in % %KE:INFO
unc_x=[0,0,0]; %noise on x; same size as x

forchange0=[];
for i=1:length(forchange)
    forchange0=[forchange0,eval(forchange{i})]; %original/predicted/supposed value
end



LB=[1 7 1037.3]'; %1,7,850 %KEINFO
UB=[2.63 20 1037.3]'; %2.63,20,1037.3 %KEINFO
perf=[];
Xopt=[];

for i=1:20
    tic
[xopt_rob fxopt M XX XXtr MMM w]=pfo_robust_v2(@optsystem_foropt_pfo,Ns,@nonlincon_pfo,LB,UB,'systematic_resampling',0,tol,uk_noise,whub_whlb,forchange,forchange0,unc_percent);
toc 
perf=[perf fxopt]
xopt_robb=[mean(XX(:,1)) mean(XX(:,2)) mean(XX(:,3))]
Xopt=[Xopt; xopt_robb] %$/h
end

end


%% FIGURE 6 (particle evolution) CIKKBE: fig6
load('fig6v1_XXtr_0515.mat')
X=XXtr;

%before u_switch
k_list=[5,20,50,70]; %KEINFO

figure
ii=1
for k=k_list
    subplot(2,2,ii)
    scatter(X{k}(:,1),X{k}(:,2),10,X{k}(:,end),'filled','o');xlabel('P_0 [bar]'),ylabel('SOR [-]'),title(sprintf('k=%2.0f',k));col=colorbar;ylabel(col,'H(\bfx\rm) [$/h]');caxis([50*10.55,200*10.55])
    ii=ii+1;
    xlim([1,2.6])
    ylim([7,10])
end

%% FIGURE 7 (counts after u _switch)
%after u_switch
k=76; %uswitch: k=76
[uu1]=unique(X{k}(:,1)); 
uu2=[];
for kk=k:length(X)
uu2=[uu2,histc(X{kk}(:,1),uu1)];
end
figure
plot(k:length(X),uu2','HandleVisibility','off'); hold on
xlabel('k [-]')
ylabel('Number of particles [-]')

uu3=find(uu2(:,end)>0); %elements where there are particles at the end
hh=plot(k:length(X),uu2(uu3,:)','LineWidth',2) %with thick line them

for ii=1:length(uu3)
    names{ii}=sprintf('\\bfx\\rm = [%2.4f,%2.4f]',X{k}(find(uu1(uu3(ii))==X{k}(:,1),1,'first'),1),X{k}(find(uu1(uu3(ii))==X{k}(:,1),1,'first'),2));
end
set(hh,{'DisplayName'},names')
lgd=legend()
lgd.Location='northwest'
lgd.Interpreter='tex'
% legend(hh,names,'Interpreter','latex')

[~,idx]=max(uu2(:,end));
% 80% of particles in:
X{end}(idx,:)
%% Result evaluation (05.09.)(CIKKBE: fig8)
%xopt: opt. with fix parameters; xopt_rob: opt. with uncertain parameters
xopt_fix=[2.03,10.03,1037.3]; %fminconból (valszeg ez a legközepe)
xopt_rob=[1.7811,9.7798,1037.3]; %xopt_robb; %robust PFO
% xopt_rob1=[1.81,9.78, 1037.3] ; %1.75
LB=[1 7 850];
UB=[2.63 20 1200]; 

%sheel1969 profit
%Parameters ($/kmol)
G_st=10.10*10.5518;
G_to=-10.36*10.5518;
G_be=-7.75*10.5518;
G_H2=0.13*10.5518;
G_meth=0.42*10.5518;
G_eth=0.68*10.5518;
G_steam1=-4.69e-3*10.5518;
G_steam2=-2.57e-5*10.5518;

E_i=[90981.4 207989.2 91515.3 103996.7 65723.3 73628.4]; %kJ/kmol
A_i=[-0.0854 13.2392 0.2961 -0.0724 -2.9344 21.2402];

% forchange={'E_i(2)','E_i(1)','A_i(2)','G_st','G_to','G_be','G_steam2','E_i(3)'};  %KE:INFO
% unc_percent=[10,10,10,15,15,15,15,0]; % in % %KE:INFO
unc_x=[0,0,0]; %noise on x; same size as x

forchange={'E_i(1)'};  %KE:INFO
unc_percent=[10]; % in % %KE:INFO

forchange0=[]; %original/predicted/supposed value
for i=1:length(forchange)
    forchange0=[forchange0,eval(forchange{i})]; %original/predicted/supposed value
end


NN=5000;
unc={};
%unc=normrnd(forchange0,unc_percent/100*forchange0,[NN,length(forchange)]); %normal
for i=1:length(forchange)
unc{i}=(rand(NN,length(forchange(i))).*2-1).*unc_percent(i)/100.*forchange0(i)+repmat(forchange0(i),length(forchange(i)),1); %uniform
end
obj_fix=[];
obj_rob=[];
obj_rob1=[];

for j=1:NN
%     E_i(2)= unc{1}(j);
%     E_i(1)= unc{2}(j);


    for i=1:length(forchange) %robust
        eval([forchange{i} ' = unc{i}(j);']); %uniform
    end

    v_i=(rand(1,3)*2-1).*unc_x/100.*(UB+LB)/2;

    obj_rob=[obj_rob;optsystem_foropt_pfo(xopt_rob+v_i)];
    obj_fix=[obj_fix;optsystem_foropt_pfo(xopt_fix+v_i)];
%     obj_rob1=[obj_rob1;optsystem_foropt_pfo(xopt_rob1+v_i)];
end

figure
hist(obj_fix,100)
xlabel('obj_{fix}')
ylabel('Histogram')
figure
hist(obj_rob,100)
xlabel('obj_{rob}')
ylabel('Histogram')

figure
hist(obj_rob-obj_fix,100)
% title('objrob-objfix')
xlabel('Objective value difference (H(\bfx\rm*_{robust})-H(\bfx\rm*_{global}))')
ylabel('Histogram of \DeltaH(\bfx\rm*)')
xline(mean(obj_rob-obj_fix),'r')
yl=ylim;
text(mean(obj_rob-obj_fix)*1.1, yl(2)*0.9,sprintf('%2.4f',mean(obj_rob-obj_fix)),'Color','r')


figure
[f,xi] =ksdensity(obj_fix);
plot(xi,f,'b')
hold on
[f,xi] =ksdensity(obj_rob);
plot(xi,f,'r')
% ksdensity(obj_rob1);
xline(mean(obj_fix),'b')
xline(mean(obj_rob),'r')
legend('Fix','Robusztus')
xlabel('Objective value')
ylabel('pdf')
%% Functions
function [c,ceq] = nonlincon_pfo(x) %=boundary(x)
global alfa_j beta_j gamma_j a_j b_j c_j L F_st0 F_be0 F_to0 UB LB
global LB UB
global H_st H_eb H_be H_tol H_steam

        %x=x.*(UB-LB)+LB;

        T_eb=800;
        P0=x(1);
        SOR=x(2);
        F_eb0=40.56/3600;
        T_steam=x(3);
        
        F_steam0=F_eb0*SOR; %kmol/s, (KE: <454 kmol/h!)

        cp_eb=alfa_j(1)+beta_j(1)*T_eb+gamma_j(1)*T_eb^2; %kJ/kmol/K, organic
        cp_steam_steam=a_j(1)+b_j(1)*T_steam+c_j(1)/T_steam^2; %kJ/kmol/K, inorganic
        cp_steam_eb=a_j(1)+b_j(1)*T_eb+c_j(1)/T_eb^2; %kJ/kmol/K, inorganic
        
        T0=(F_eb0*cp_eb*T_eb+50.4/3600*cp_steam_eb*T_eb+(F_steam0-50.4/3600)*cp_steam_steam*T_steam)/(F_eb0*cp_eb+(F_steam0-50.4/3600)*cp_steam_steam+50.4/3600*cp_steam_eb); %kJ/s / kJ/s/K = K

        c(1)=F_steam0-454/3600;
        c(2)=850-T0;
        c(3)=T0-925;
        ceq = [];
        %c=[];
        
end

function obj=optsystem_foropt_pfo(x) %=objective(x)
global ment ment2 ment6
global alfa_j beta_j gamma_j a_j b_j c_j L  
global LB UB ment2 ment6
global H_st H_eb H_be H_tol H_steam


%x=x.*(UB-LB)+LB;

T_eb=800;
P0=x(1);
SOR=x(2);
F_eb0=40.56/3600;
T_steam=x(3);

F_st0=0.67/3600*(F_eb0*3600/36.87); %kmol/s, styrene
F_be0=0.11/3600*(F_eb0*3600/36.87); %kmol/s, benzene
F_to0=0.88/3600*(F_eb0*3600/36.87); %kmol/s, toluene

F_steam0=F_eb0*SOR; %kmol/s, (KE: <454 kmol/h!)

cp_eb=alfa_j(1)+beta_j(1)*T_eb+gamma_j(1)*T_eb^2; %kJ/kmol/K, organic
cp_steam_steam=a_j(1)+b_j(1)*T_steam+c_j(1)/T_steam^2; %kJ/kmol/K, inorganic
cp_steam_eb=a_j(1)+b_j(1)*T_eb+c_j(1)/T_eb^2; %kJ/kmol/K, inorganic

T0=(F_eb0*cp_eb*T_eb+50.4/3600*cp_steam_eb*T_eb+(F_steam0-50.4/3600)*cp_steam_steam*T_steam)/(F_eb0*cp_eb+(F_steam0-50.4/3600)*cp_steam_steam+50.4/3600*cp_steam_eb); %kJ/s / kJ/s/K = K

ment2=[ment2;T0]; %saving the inlet T of the eb+steam mixture, (KE: <925 K)
ment6=[ment6;F_steam0*3600];

[z x]=ode23s(@systemmodel,0:0.1:L,[F_eb0 F_st0 F_be0 F_to0 0 0 F_steam0 0 0 0 T0 P0]');
x(:,1:10)=x(:,1:10)*3600; %kmol/s --> kmol/h
%ment=[ment;x(end,1:10)]; %saving the F_j mole flows

obj1=obj_yee(x);
% obj2=obj_clough(x,T_steam);
obj3=obj_sheel(x,T_steam);
obj1(2:3)=obj1(2:3)*100;
obj=obj3; %obj1(4); %PFO: maximization
end


function objval=obj_yee(x) %yee2003
global H_st H_eb H_be H_tol H_steam ment7
    profit=x(end,2)*H_st+x(end,3)*H_be+x(end,4)*H_tol-(x(1,1)-x(end,1))*H_eb-x(end,7)*H_steam; %$/h, profit
    S_st=(x(end,2)-x(1,2))/(x(1,1)-x(end,1)); %selestivity
    F_st=x(end,2); %amount of styrene produced
    Y_st=(x(end,2)-x(1,2))/x(1,1); %yield of styrene

    ment7=[ment7;x(end,2) x(end,3) x(end,4) x(1,1)-x(end,1) x(end,7)];
    objval=[F_st,S_st,Y_st,profit];
end

function objval=obj_sheel(x,T_steam)
global G_st G_to G_be G_H2 G_meth G_eth G_steam1 G_steam2
profit=x(end,2)*G_st+x(end,4)*G_to+x(end,3)*G_be+x(end,8)*G_H2+x(end,6)*G_meth+x(end,5)*G_eth+(G_steam1+G_steam2*T_steam)*x(1,7); %$/h
objval=profit;
end


function dx_dz = systemmodel(z,x)
    global G0  Nc Nr L rho_b D_p e A_t nu_ij E_i A_i M_j G a_i b_i alfa_j beta_j gamma_j a_j b_j c_j ment ment1 F_st0 F_be0 F_to0
    %x=[Fj (j=1,...,Nc) [kmol/s],T [K],P [bar]]
    %j=[(eb),(st),(be),(tol),(eth),(meth),steam,H2,CO,CO2] 

    Keb=exp(-(122725-126.3*x(Nc+1)-0.002194*x(Nc+1)^2)/(8.314*x(Nc+1))); %bar
    %x(1:Nc)=abs(x(1:Nc));
    p_j=x(1:Nc)/sum(x(1:Nc))*x(Nc+2); %bar
    %p_j(p_j<0)=0;
    k_i=exp(A_i-E_i/8.314/x(Nc+1)); %kmol/kg cat./s/bar^n

    r_i(1)=k_i(1)*(p_j(1)-p_j(2)*p_j(8)/Keb);
    r_i(2)=k_i(2)*p_j(1);
    r_i(3)=k_i(3)*p_j(1)*p_j(8);
    r_i(4)=k_i(4)*p_j(7)*p_j(5)^(0.5);
    r_i(5)=k_i(5)*p_j(7)*p_j(6);
    r_i(6)=k_i(6)*x(Nc+2)/x(Nc+1)^3*p_j(7)*p_j(9); %kmol/s/kg cat.

   % ment=[ment; r_i(6)];
  
    R_j=r_i*nu_ij; %kmol/kg cat./s  %sor(1x10)=sor(1x6)xoszlop(6x10)
    dx_dz(1:Nc)=rho_b*A_t*R_j; %kmol/s/m

    %  ment1=[ment1; R_j(5:6)];

    dH_i=a_i+b_i*x(Nc+1); %kJ/kmol
    cp_j(1:6)=alfa_j+beta_j*x(Nc+1)+gamma_j*x(Nc+1)^2; %kJ/kmol/K, organic
    cp_j(7:10)=a_j+b_j*x(Nc+1)+c_j/x(Nc+1)^2; %kJ/kmol/K, inorganic

    dx_dz(Nc+1)=sum(-dH_i.*r_i*rho_b*A_t,2)/sum(x(1:Nc)'.*cp_j,2); %(kJ/s/m) / (kJ/s/K) = K/m

    M_av=sum(x(1:Nc)'.*M_j)/sum(x(1:Nc)); %kg/kmol
    rho_G=10^5*x(Nc+2)*M_av/x(Nc+1)/(8.314*1000); %kg/m^3 ([R]=J/kmol/K)
    mu_G=3e-5; %Pas %1.6e-5 
    G=sum(x(1:Nc)'.*M_j); %kg/s
    G0=G/A_t; %kg/m^2/s

    dx_dz(Nc+2)=-10^(-5)*(1-e)*G0/D_p/e^3/rho_G*(150*(1-e)*mu_G/D_p+1.75*G0); %10^(-5)*Pa=bar
    dx_dz=dx_dz';
end

function obj=optsystem_foropt_fm(x)
global ment  ment2 ment6
global alfa_j beta_j gamma_j a_j b_j c_j L F_st0 F_be0 F_to0 F_eb0 cp_eb T_eb
global LB UB
global H_st H_eb H_be H_tol H_steam C_st C_be C_to C_h1 C_h2
global G_st G_to G_be G_H2 G_meth G_eth G_steam1 G_steam2

x=x.*(UB-LB)+LB;

T_eb=800;
P0=x(1);
SOR=x(2);
F_eb0=40.56/3600;
T_steam=x(3);

F_st0=0.67/3600*(F_eb0*3600/36.87); %kmol/s, styrene
F_be0=0.11/3600*(F_eb0*3600/36.87); %kmol/s, benzene
F_to0=0.88/3600*(F_eb0*3600/36.87); %kmol/s, toluene

F_steam0=F_eb0*SOR; %kmol/s, (KE: <454 kmol/h!)

cp_eb=alfa_j(1)+beta_j(1)*T_eb+gamma_j(1)*T_eb^2; %kJ/kmol/K, organic
cp_steam_steam=a_j(1)+b_j(1)*T_steam+c_j(1)/T_steam^2; %kJ/kmol/K, inorganic
cp_steam_eb=a_j(1)+b_j(1)*T_eb+c_j(1)/T_eb^2; %kJ/kmol/K, inorganic

T0=(F_eb0*cp_eb*T_eb+50.4/3600*cp_steam_eb*T_eb+(F_steam0-50.4/3600)*cp_steam_steam*T_steam)/(F_eb0*cp_eb+(F_steam0-50.4/3600)*cp_steam_steam+50.4/3600*cp_steam_eb); %kJ/s / kJ/s/K = K

[z x]=ode23s(@systemmodel,0:0.1:L,[F_eb0 F_st0 F_be0 F_to0 0 0 F_steam0 0 0 0 T0 P0]');
x(:,1:10)=x(:,1:10)*3600; %kmol/s --> kmol/h
%ment=[ment;x(end,1:10)]; %saving the F_j mole flows

obj1=obj_yee(x);
% obj2=obj_clough(x,T_steam);
obj3=obj_sheel(x,T_steam);
% obj4=obj_utility(x,T_steam);
obj1(2:3)=obj1(2:3)*100;
obj=-obj3; %KE: INNEN
end

function [c,ceq] = nonlincon_fm(x)
global alfa_j beta_j gamma_j a_j b_j c_j L  UB LB
global LB UB
global H_st H_eb H_be H_tol H_steam
global G_st G_to G_be G_H2 G_meth G_eth G_steam1 G_steam2

        x=x.*(UB-LB)+LB;

        T_eb=800;
        P0=x(1);
        SOR=x(2);
        F_eb0=40.56/3600;
        T_steam=x(3);
        
        F_steam0=F_eb0*SOR; %kmol/s, (KE: <454 kmol/h!)

        cp_eb=alfa_j(1)+beta_j(1)*T_eb+gamma_j(1)*T_eb^2; %kJ/kmol/K, organic
        cp_steam_steam=a_j(1)+b_j(1)*T_steam+c_j(1)/T_steam^2; %kJ/kmol/K, inorganic
        cp_steam_eb=a_j(1)+b_j(1)*T_eb+c_j(1)/T_eb^2; %kJ/kmol/K, inorganic
        
        T0=(F_eb0*cp_eb*T_eb+50.4/3600*cp_steam_eb*T_eb+(F_steam0-50.4/3600)*cp_steam_steam*T_steam)/(F_eb0*cp_eb+(F_steam0-50.4/3600)*cp_steam_steam+50.4/3600*cp_steam_eb); %kJ/s / kJ/s/K = K


        c(1)=F_steam0-454/3600;
        c(2)=850-T0;
        c(3)=T0-925;
        ceq = [];
end

