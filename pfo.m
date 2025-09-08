function [xopt fxopt M XX XXtr MMM w]=pfo(objective,Ns,boundary,LB,UB,resampling_strategy,fi_exp,tol,uk_noise,weightrate)
global ment2 ment6 bbb whub_whlb
%Particle Filter Optimisation algorithm. Resampling with replacement is
%used, and the particles are forced to move by the transition kernel, thus
%ensure the exploration of the decision variable space. The system noise
%(uk)
%was handled by a trajectory like (UB-LB)/(a*k^b) and decrease in time. The
%distribution of the 'vk' observation noise (fi) can be used as uniform
%or exponential, too. It can be set by the fi_exp switch.
%The arguments of the function:
%objective: the objective function which inputs are the decision variables
%and its output (only one!) is the objective value.
%Ns:number of the particles.
%LB,UB: lower and upper bounds of the decision variables. Give them in a column vector
%resampling_strategy: it can be 'multinomial_resampling' or
%'systematic_resampling'.
%fi_exp: it can have two value: 0 if uniformly distributed fi wanted to
%used and 1, if exponentially distributed.
%tol: the tolerance of the stop criteria.
%uk_noise: it has to be a 1x2 vector which determines the coefficients of
%uk like [a b].
% 
% Nx: number of decision variables

% tol=0.0001; %KE: hangolandó
% uk_noise=[65/30 1.25]; %(Hub-Hlb)/noise, KE: hangolandó,50/30
% whub_whlb=300; %w(Hub)/w(Hlb); KE: hangolandó
 whub_whlb= weightrate;

Nx=size(LB,1);
k=1;
ment2=[];
ment6=[];
stop=0;
clear M X xopt uk fx y w X_res w_res MM MMM


if fi_exp==1 %exponential
    fi=@fi3;
    vk=@vk3;
else %uniform
    fi=@fi;
    vk=@vk;
end
szaml=0;

while stop==0 %stopping cr.
    i=1;
    while i<=Ns
         if k==1
            %Initialization
            X{k}(i,:)=rand(1,Nx).*(UB-LB)'+LB'; % Ns x Nx
         else
            %Importance sampling (transition)
%                 if k==2
%                 uk(i,:)=mvnrnd(0,1,Nx)'.*(UB-LB)'./(uk_noise(1)*k^uk_noise(2));
%                 else
                %uk(i,:)=randn(1,Nx).*(UB-LB)'./(100);
uk(i,:)=mvnrnd(zeros(Nx,1),eye(Nx).*((UB-LB)'./(uk_noise(1)*k^uk_noise(2)))')';    %eye, multi , uk_noise=[1.2 3.25];

    %       uk(i,:)=mvnrnd(0,1,Nx)'.*(UB-LB)'./(uk_noise(1)*k^uk_noise(2));  %exponential I.  , uk_noise=[1.5 1.65]; 
  %  uk(i,:)=mvnrnd(0,1,Nx)'.*(UB-LB)'./(uk_noise(1)*uk_noise(2)^k);
%     %exponential II.

                
                % %                 if k-1>2
% %                 uk(i,:)=mvnrnd(0,1,Nx)'.*(UB-LB)'./(uk_noise(3)*(mean(M(k-1)))^uk_noise(4));
% %                 end
%                 end
            X{k}(i,:)=X_res{k-1}(i,:)+uk(i,:); %with transition
            %X{k}=X_res{k-1}; %without transition
         end
         % Checking the boundaries
         %if  prod(boundary(X{k}(i,:))<0)==0 || prod(X{k}(i,:)<UB)==0 || prod(X{k}(i,:)>LB)==0 %if one condition is false, repeat the sampling
         if  any([boundary(X{k}(i,:))>0,X{k}(i,:)>UB',X{k}(i,:)<LB']) %if one condition is false, repeat the sampling
         continue
         end
         %Observation construction
         fx{k}(i)=objective(X{k}(i,:)); %KE: zhou2011 szerint Hx, nálam fx a jelölés
         i=i+1;    
    end
     if max(fx{k})==min(fx{k})
        stop(k)=1
        continue
     end

    y(k)=max(fx{k})-vk(max(fx{k}),min(fx{k}));

   

    if k>1
    if y(k)<y(k-1)
        y(k)=y(k-1);
    end
    end

    % Bayes' updating (generating the new normalized weights)
    if k==1
        sum_wk=sum(fi(fx{k}-repmat(y(k),size(fx{k},1),size(fx{k},2)),max(fx{k}),min(fx{k})));
        w{k}=fi(fx{k}-repmat(y(k),size(fx{k},1),size(fx{k},2)),max(fx{k}),min(fx{k}))/sum_wk; %1xNs; unnormalized
    else
        sum_wk=sum(fi(fx{k}-repmat(y(k),size(fx{k},1),size(fx{k},2)),max(fx{k}),min(fx{k})).*w_res{k-1});
        if sum_wk==0; %KE: arra az esetre, ha az új yk kisebb mint az előző, és annál jobb részecske nincs már
            szaml=szaml+1;
            if max(fx{k})-min(fx{k})<max(fx{k})/100
                stop(k)=1;
                continue
            else 
                if szaml==10
                    stop(k)=1;
                    fx{k}=fx{k-1};
                    X{k}=X{k-1};
                end
                continue
            end
        else 
            szaml=0;
        end
        w{k}=fi(fx{k}-repmat(y(k),size(fx{k},1),size(fx{k},2)),max(fx{k}),min(fx{k})).*w_res{k-1}/sum_wk; %1xNs; unnormalized
    end
    w{k} =w{k}./sum(w{k});
    %Resampling
        %resampling with replacement
        [X_res{k}, w_res{k}, idx{k}] = resample(X{k}', w{k}', resampling_strategy); %x row, w column 
        X_res{k}=X_res{k}';
     %Define the stop criteria
     M(k)=0;
     if k>10
         %M(k)=mean(X{k});
         %M(k)=mean(fx{k});

        M(k)=sum(fx{k}.*w{k})/sum(w{k}); %k>2!!!
         stop(k)=any([all([abs(M(k)-M(k-1))<tol,uk_noise(1)*k^uk_noise(2)>0]),max(fx{k})-min(fx{k})<1]);%*ones(1,Nx);
       %disp(M(k));

%          %liu2016_2 (k>10!!!)
%         MM(k)=max(fx{end});  
%         stop(k)=~ismember(max(MM),MM(k-9:k));

         %stop(k)=sum(std(X_res{k}))>0; 
         
     else
         stop(k)=0;
     end
MMM(k)=max(fx{k});
 k=k+1;
end

%Case1: xopt is the expected value of the particles
% xopt=sum(X{end}.*repmat(w{end}',1,size(X{end},2)))./sum(w{end});
%KE: így hogy benn van az intervallumos leállási feltétel continue-val msotmár nem jó
%ez, emrt az utolsó k-nak (lehet?) nem lesznek súlyai

%Case2: xopt is the x that belongs to the maxima of fx{end}
% [~,id]=max(fx{end});
% xopt=X{end}(id,:);


%Case3: xopt is the maxima of fx from all iterations
MM=MMM %HA saját, HA lui, törlés! (van MM)
[~,id2]=max(MM); %liu2016_2 %id2: # of the run
[~,id3]=max(fx{id2}); %id3: # of the particle
xopt=X{id2}(id3,:);

% eva = evalclusters(X{end},'kmeans','CalinskiHarabasz','klist',1:10)
% eva = evalclusters(X{end},'kmeans','gap','klist',1:10)
% %finding the local optimas with clustering
% [classes,xopt]=kmeans(X{end},eva.OptimalK)


fxopt=objective(xopt);
XX=[X{end} fx{end}'];
for i=1:size(fx,2)
XXtr{i}=[X{i} fx{i}'];
end
if exist('MM')
M=MM;
else
%     M=zeros(1,k);
end
end





function vk_rand=vk(Hub,Hlb) %generate a random vk according to fi
    pd=makedist('Uniform','Lower',0,'Upper',Hub-Hlb);
    vk_rand=random(pd);
end

function pval=fi(Hx_yk,Hub,Hlb) %evaluate fi() in Hx-yk
    pd=makedist('Uniform','Lower',0,'Upper',Hub-Hlb);
    pval=pdf(pd,Hx_yk);
end

function vk_rand=vk3(Hub,Hlb) %generate a random vk according to fi (exponential distribution)
    global whub_whlb
    
    pd=makedist('Exponential','mu',(Hub-Hlb)/log(whub_whlb)); %if log(1) it is the same as uniform
    pd=truncate(pd,0,Hub-Hlb);
    vk_rand=Hub-Hlb-random(pd); %nagy vk; yk<0;
    
   
end

function pval=fi3(Hx_yk,Hub,Hlb) %evaluate fi() in Hx-yk
    global whub_whlb
    pd=makedist('Exponential','mu',(Hub-Hlb)/log(whub_whlb));
    pd=truncate(pd,0,Hub-Hlb);
    pval=pdf(pd,Hub-Hlb-Hx_yk);
end

function vk_rand=vk2(Hub,Hlb) %generate a random vk according to fi (exponential distribution)
    global whub_whlb
    pd=makedist('Exponential','mu',(Hub-Hlb)/log(whub_whlb)); %if log(1) it is the same as uniform
    vk_rand=random(pd);
end

function pval=fi2(Hx_yk,Hub,Hlb) %evaluate fi() in Hx-yk
    global whub_whlb
    pd=makedist('Exponential','mu',(Hub-Hlb)/log(whub_whlb));
    pval=pdf(pd,Hx_yk);
end


