function main_ecoevo(varargin)

%This is the main function for the eco-evolutionary dynamics analysis
%The function takes input arguments for: 
% -  Prey abundance (N0) in habitat 1,2, and 3 (total 3 input arguments, arguments one, two, and three)
% - Prey dispersal (mN, argument four)
% - Prey niche width (sigma_a, argument five)
% - Predator abundance (P0) in habitat 1, 2, and 3 (total 3 input arguments, arguments six, seven and eight)
% - Predator dispersal (mP, argument nine). 
% - Predator niche width (sigma_b, argument ten)
% - Index for replicate (argument eleven, will show in saved output)

%Note if no arguments are assigned then default values (see below) will be
%run. 

%Note if input arguments are assigned all 11 arguments listed above needs to be assigned
%or else default values will be used.

rand('seed',fix(sum(1000*clock)));
randn('seed',fix(sum(2000*clock)+1));

name=['mainout.mat'];


%% Set up the default paramters and innitial conditions of the ecological model

%------Set habitat variables-----------
%Parameters for the resource distributions in habitats A and B
K0=[10000; 12000; 13000]; %Max K in habitats (here one value i.e. same K in both habitats) 
U=[0 1 2]; %Resource peeks in the two habitats
sigma_K=[1 1 1]; %Width of the resource distribution (here one value i.e. same in both habitats)

%-----Set competitor variables-------
N0=[2 0 0;
    0 2 0;
    0 0 2]; %Row i in matrix denote abundance for species i in habitat j (collumns) 
 
V0=[0 1 2]; %Trait, each species initally
V=V0; %Trait, each species
prey_sp_id =zeros(size(V)); %Vector for species id
prey_sp_id_counter = 0; 
r=1; %intrisic growth rate
mN=0.05; %Consumer dispersal
sigma_a=0.5; %Width of the consumption kernel (niche width)
mut_N=0.01; %mutation rate
sigma_mut_N = 0.02; % standard deviation of consumer mutations

%-----Set predator variables ---------
%Set predator variables
P0=[1 0 0;
    0 1 0;
    0 0 1]; %Row i in matrix denote abundance for species i in habitat j (collumns) 

Z0=[0 1 2]; %Trait, each element, one species
Z = Z0;
pred_sp_id =zeros(size(Z)); %Vector for species id
pred_sp_id_counter = 0;
d=-0.2; %intrisic death rate
mP=0.05; %Predator dispersal
sigma_b=0.4; %niche width of the predator
mut_P=0.1; %mutation rate
sigma_mut_P = 0.03; % standard deviation of predator mutations
bmax=0.0001; %Max attack rate (i.e. attack rate when complete trait match between pred and prey)
cP=0.3; %Conversion from competitor to predator

%---------Implementation parameters-------
epsi=1; %Extinction threshold
t_evo=1; %Evolutonary time, first  

%%
%If input arguments were assigned 

if nargin == 11
    
    
    disp('User input for N0:')
    N0=[varargin{1} 0 0;
        0 varargin{2} 0;
        0 0 varargin{3}]
    disp('User input for mN:')
    mN = varargin{4}
    disp('User input for sigma_a:')
    sigma_a= varargin{5}
    
    disp('User input for P0:')
    P0=[varargin{6} 0 0;
        0 varargin{7} 0;
        0 0 varargin{8}]
    disp('User input for mP:')
    mP = varargin{9}
    disp('User input for sigma_b:')
    sigma_b= varargin{10}
    
    disp('User input for replicate index:')
    rep_ind=varargin{11}
    disp('Output name of this run:')
    name=['mainout' num2str(rep_ind) '.mat']

    
    
    disp('This program will continute in 20 seconds...')
    pause(20)
    
else 
    disp('Note: Eleven input arguments must be assigned or default parameter values will run')
    disp('This program will continue with default values in 10 seconds...')
    
    pause(10)
    
end 


%%
%Rearange N0 and P0 such that it fits the ODE system
[m_cons n]=size(N0);
[m_pred dumy]= size(P0);

tmp=[N0;P0];
[m_tot dumy]=size(tmp);

abun_vec = reshape(tmp,1,m_tot*n)';

%%
%Save the parameters in a struct that will be piped into functions later
parameters.abun_vec = abun_vec;
parameters.V=V0; %Initial conditions traits
parameters.Z=Z0; %Initial conditions traits
parameters.r=r; %Consumer growth
parameters.d=d; %Predator death
parameters.K0=K0; %Max K
parameters.bmax = bmax; %Max predator consumption
parameters.sigma_a=sigma_a; %Width consumption kernel
parameters.sigma_b = sigma_b; %Width predation kernel
parameters.mN=mN; %Consumer dispersal
parameters.mP=mP; %Predator dispersal
parameters.cP=cP; %Predator conversion rate
parameters.U=U; %Resource peaks
parameters.sigma_K=sigma_K; %Resource width
parameters.m_cons = m_cons; %No. consumers
parameters.m_pred = m_pred; %No. consumers
parameters.n = n; %No. habitats
parameters.m_tot =m_tot; %Total number of species

%%
%Plot the resource landscape and species niches for the initial conditions
PlotResLand(parameters)

%%
%Compute and plot the populaton dynamics
[t,y,K_a,alpha,a]=pop_dyn_func(parameters,100000);
figure(2)
subplot(2,1,1)
plot( t , y)
% legend('N_1 in A','N_2 in A','N_3 in A','N_1 in B','N_2 in B','N_3 in B','N_1 in C','N_2 in C','N_3 in C')
title('Equilibrium and dynamics of the seeded community')
axis square

%Plot species dynamics
PlotSpecDyn(y,t,abun_vec,parameters)

%%
%Compute initial fitness landscape
mut_range=-5:0.01:5; %a range of trait values that will be used to compute fitness lanscape
mut_hab = 1:n; %Habitats for which we test mutation fitness
abun_vec=y(end,:)'; %Collect equilibrium population sizes 
N_tot = reshape(abun_vec,m_tot,n); %Reshape abun_vec to abundanse matrix where rows are species and collumns are habitats 
N=N_tot(1:m_cons,:);
P=N_tot(m_cons+1:end,:);

% PlotInitFitLand(mut_range,mut_hab,abun_vec,N_tot,N,P,V,Z,parameters)

%%
%Before we go into the evolutionary loop we do three things:

%First remove populations that may have gone extinct at the first
%calculation of the equilibrium

%Second we check if the initial community is constituted of more than one
%species. This is done by checking for gaps in the trait distributions. 

%Third we save the initial community trait distribution, 
%abundance and species Id at equilibrium

%Perge the system from extinct populations
[N V P Z m_cons n m_pred N_tot m_tot abun_vec prey_sp_id pred_sp_id]=PergeFunc(N,P,V,Z,prey_sp_id,pred_sp_id);

%Save first set of information for the output
prey_dist_data{1,1}=1; %time
prey_dist_data{1,2}=[V; prey_sp_id]; %matrix for trait distribution and species id
prey_dist_data{1,3}=N; %Abundance for each morph in each habitat

prey_phylo_data(1,1)=0; %species id
prey_phylo_data(1,2)=0; %origin
prey_phylo_data(1,3)=0; %time first registered

pred_dist_data{1,1}=1; %time
pred_dist_data{1,2}=[Z; pred_sp_id]; %matrix for trait distribution and species id
pred_dist_data{1,3}=P; %Abundance for each morph in each habitat

pred_phylo_data(1,1)=0; %species id
pred_phylo_data(1,2)=0; %origin
pred_phylo_data(1,3)=0; %time first registered

com_info_data{1,1}=1; %time
com_info_data{1,2}=K_a; %K_A
com_info_data{1,3}=alpha; %alpha
com_info_data{1,4}=a; %a

%Ceck for separate species in the starting conditions
check_PreyPhylo_flag=1;
check_PredPhylo_flag=1;
for i=1:length(V) %Do this multiple times to make sure all gaps are considered
    [prey_sp_id, pred_sp_id, prey_sp_id_counter, pred_sp_id_counter, prey_phylo_data, pred_phylo_data]=PhyloCheck(V, prey_sp_id, Z, pred_sp_id,prey_sp_id_counter,pred_sp_id_counter, sigma_mut_N, sigma_mut_P, check_PreyPhylo_flag,check_PredPhylo_flag,t_evo, prey_phylo_data, pred_phylo_data);
end



%%
%HERE WE START THE EVOLUTOINARY TIMELINE, LOOPING OVER EVOLUTIONARY TIME
%Compute the rate w for the emergence of mutants, pick mutating population and create mutant
%This rate is dependent on intrisic birth rates(see Ito and Dieckmann)
%We have the birth (for consumers) and deaths (predators) in our
%growth_death_vec. As births=deaths at equilibrium I simply take the
%absolute of deaths in this vector to get birth of predators. 
tmp1 = N.*r*mut_N;
tmp2 = P.*abs(d)*mut_P;
tmp3 = [tmp1;tmp2];
w = reshape(tmp3,1,m_tot*n)';
w_tot=sum(w);


%Loop over evolutionary time
count=0;
while t_evo < 6100
    count=count+1;
    check_PreyPhylo_flag=0;
    check_PredPhylo_flag=0;
    
    %Pick the population that should mutate according to w/w_tot
    mut_pop=[];
    tmp=1:length(w);
    mut_pop=tmp(find(rand<cumsum(w/w_tot),1,'first'));
    [mut_pop mutHab] = find(N_tot==abun_vec(mut_pop)); %Get mutpop and mutHab
    %More than one population can be picked so we need to randomly pick
    %only one, here the first one
    mut_pop=mut_pop(1);
    mutHab=mutHab(1);
        
    %Find out if the mutating population is a consumer or a predator
    %This can be done by looking at the groth_death_vector
    
    if mut_pop<=m_cons %if the mutating species is a competitor

        disp('Consumer mutated')
        
        %Get the mutant trait and species id
        mut_V=V(mut_pop) + sigma_mut_N*randn;
        mut_spId=prey_sp_id(mut_pop);
     
        %Update evolutionary time according to ITO & DIECKMANN 
        t_evo=t_evo+1; %t_evo+(-(1/w_tot)*log(rand)); 
        
        %Set up mutant
        N_tmp = N; 
        N_tmp(end+1,:)=[0,0,0]; %Ad a mutant 
        N_tmp(end,mutHab)=1; %Introduce one mutant individual in the correct habitat
        N_tmp(mut_pop,mutHab) = N_tmp(mut_pop,mutHab)-1; %Remove a individual from the mutating species
        V_tmp = V; %Set up trait vector that includes mutant
        V_tmp(end+1)=mut_V;
        
        [m_cons_tmp n]=size(N_tmp); %no. rows represent no of species and no collumns represent no habitats
        N_tot_tmp=[N_tmp;P];
        abun_vec_tmp = reshape(N_tot_tmp,1,(m_cons_tmp+m_pred)*n)';

        parameters.N_tot = N_tot_tmp;
        parameters.N = N_tmp;
        parameters.P = P;
        parameters.abun_vec=abun_vec_tmp; %Collect equilibrium population sizes 
        parameters.mutHab=mutHab;
        parameters.V=V_tmp;
        parameters.Z = Z; 

        %Compute invasion fitness            
        [t,y,mutFit_loc_cons mutFit_glob_cons mutFit_sim_cons mutInvFit_cons]=cons_fitness_func(parameters,100);            

        %Invade with probability that is based on fitness
        tmp=rand; %draw a random number between 0-1

        if mutFit_sim_cons>0; %Allow all mutants with positive invasion fitness to invade 
        %if tmp < mutFit_sim_cons/r; %Let mutant invade with some fitness dependent probability
            
            disp('Consumer mutant invaded')
            check_PreyPhylo_flag =1;
            %Save N_tmp, abun_vec_tmp and v_tmp. These will be used to compute
            %equilibrium if mutual invasibility exist
            N_mutinv = N_tmp; abun_vec_mutinv=abun_vec_tmp; V_mutinv = V_tmp; 
            
            %Replace the resident with the mutant, compute new equilibrium and 
            %test for mutual invasibility
            parameters.abun_vec=abun_vec; %Collect original equilibrium population sizes, not the tmp one from the fitnes calculation 
            V_tmp = V; %Get original V-vector but switch to mutant trait, see next line 
            V_tmp(mut_pop)=mut_V;
            parameters.V=V_tmp;
                        
            [t,y,K_a,alpha,a]=pop_dyn_func(parameters,10000);

            %Save data for next evolutionary step
            V_4next = V_tmp;
            prey_sp_id4next=prey_sp_id;
            y(y<0) = 0;%Remove negative abundance values, even though they should not occur according to the ODE solver settings
            y_4next = y; 
            t_4next = t;
            
            %Compute mutual invasibility
            %Set up resident as a mutant
            abun_vec_tmp =y(end,:)'; %Collect equilibrium population sizes 
            N_tot_tmp = reshape(abun_vec_tmp,m_tot,n); %Reshape abun_vec to abundanse matrix where rows are species and collumns are habitats 
            N_tmp=N_tot_tmp(1:m_cons,:);
            N_tmp(end+1,:)=[0,0,0]; %Ad a mutant (here the resident) 
            N_tmp(end,mutHab)=1; %Introduce one mutant individual in the correct habitat
            N_tmp(mut_pop,mutHab) = N_tmp(mut_pop,mutHab)-1; 
            V_tmp(end+1)=V(mut_pop);
            
            P_tmp = N_tot_tmp(m_cons+1:end,:);
            
            [m_cons_tmp n]=size(N_tmp); %no. rows represent no of species and no collumns represent no habitats
            N_tot_tmp=[N_tmp;P_tmp];
            abun_vec_tmp = reshape(N_tot_tmp,1,(m_cons_tmp+m_pred)*n)';

            parameters.N_tot = N_tot_tmp;
            parameters.N = N_tmp;
            parameters.P = P_tmp;
            parameters.abun_vec=abun_vec_tmp; %Collect equilibrium population sizes 
            parameters.mutHab=mutHab;
            parameters.V=V_tmp;

            %Compute resident invasion fitness            
            [t,y,mutFit_loc_cons mutFit_glob_cons mutFit_sim_cons mutInvFit_cons]=cons_fitness_func(parameters,100);
                        
             if mutFit_sim_cons > 0 %If mutual invasibility exist
               
                disp('Mutual invasibility exist')
                %Compute new equilibrium with both resident and mutant present
                parameters.N_tot = [N_mutinv; P];
                parameters.N = N_mutinv;
                parameters.P = P;
                parameters.abun_vec=abun_vec_mutinv; %Collect equilibrium population sizes 
                parameters.V=V_mutinv;
                
                [t,y,K_a,alpha,a]=pop_dyn_func(parameters,10000);

                %Update data for next evolutionary step
                 %Save data for next evolutionary step
                V_4next = V_mutinv;
                prey_sp_id4next=[prey_sp_id prey_sp_id(mut_pop)];
                y(y<0) = 0;%Remove negative values, even though they should not occur according to the ODE solver settings
                y_4next = y; 
                t_4next = t;
             end
             
            %Update abundance and trait data
            V = V_4next;
            prey_sp_id=prey_sp_id4next;
            m_cons = length(V);
            abun_vec=y_4next(end,:)'; %Collect equilibrium population sizes 
            m_tot=m_cons+m_pred;
            N_tot = reshape(abun_vec,m_tot,n); %Reshape abun_vec to abundanse matrix where rows are species and collumns are habitats 
            N=N_tot(1:m_cons,:);
            P=N_tot(m_cons+1:end,:);
            
           
            
            %Perge the system from populations with abudnance < 1
            [N V P Z m_cons n m_pred N_tot m_tot abun_vec prey_sp_id pred_sp_id]=PergeFunc(N,P,V,Z,prey_sp_id,pred_sp_id);

            
            %Recompute the rate w 
            tmp1 = N.*r*mut_N;
            tmp2 = P.*abs(d)*mut_P;
            tmp3 = [tmp1;tmp2];
            w = reshape(tmp3,1,m_tot*n)';
            w_tot=sum(w);
            
                        
        end % end of if invasion fitness is positive test
        
        
    elseif mut_pop>m_cons %if the mutating species is a predator and the predator exist
        
        disp('Predator mutated')
        mut_pop = mut_pop-m_cons;
        
        %Get the mutant trait and species id
        mut_Z=Z(mut_pop) + sigma_mut_P*randn;
        mut_spId = pred_sp_id(mut_pop);
       
        %Update evolutionary time according to ITO & DIECKMANN 
        t_evo=t_evo+1;%t_evo+(-(1/w_tot)*log(rand)); 
        
        %Set up mutant
        P_tmp = P; 
        P_tmp(end+1,:)=[0,0,0]; %Ad a mutant 
        P_tmp(end,mutHab)=1; %Introduce one mutant individual in the correct habitat
        P_tmp(mut_pop,mutHab) = P_tmp(mut_pop,mutHab)-1; %Remove a individual from the mutating species
        Z_tmp = Z; %Set up trait vector that includes mutant
        Z_tmp(end+1)=mut_Z;
        
        [m_pred_tmp n]=size(P_tmp); %no. rows represent no of species and no collumns represent no habitats
        N_tot_tmp=[N;P_tmp];
        abun_vec_tmp = reshape(N_tot_tmp,1,(m_cons+m_pred_tmp)*n)';

        parameters.N_tot = N_tot_tmp;
        parameters.N = N;
        parameters.P = P_tmp;
        parameters.abun_vec=abun_vec_tmp; %Collect equilibrium population sizes 
        parameters.mutHab=mutHab;
        parameters.Z=Z_tmp;
        parameters.V = V; 

        %Compute invasion fitness            
        [t,y,mutFit_loc_pred mutFit_glob_pred mutFit_sim_pred mutInvFit_pred]=pred_fitness_func(parameters,100);            

        %Invade with probability that is based on fitness
        tmp=rand; %draw a random number between 0-1

        if mutFit_sim_pred>0; %Allow all mutants with positive invasion fitness to invade 
        %if tmp < mutFit_sim_pred/r; %Let mutant invade with some fitness dependent 

            disp('Predator mutant invaded')
            check_PredPhylo_flag=1;
            %Save P_tmp, abun_vec_tmp and Z_tmp. These will be used to compute
            %equilibrium if mutual invasibility exist
            P_mutinv = P_tmp; abun_vec_mutinv=abun_vec_tmp; Z_mutinv = Z_tmp; 
            
            %Replace the resident with the mutant, compute new equilibrium and 
            %test for mutual invasibility
            parameters.abun_vec=abun_vec; %Collect original equilibrium population sizes, not the tmp one from the fitnes calculation 
            Z_tmp = Z; %Get original V-vector but switch to mutant trait, see next line 
            Z_tmp(mut_pop)=mut_Z;
            parameters.Z=Z_tmp;
                        
            [t,y,K_a,alpha,a]=pop_dyn_func(parameters,10000);

            %Save data for next evolutionary step
            Z_4next = Z_tmp;
            pred_sp_id4next=pred_sp_id;
            y(y<0) = 0;%Remove negative abundance values, even though they should not occur according to the ODE solver settings
            y_4next = y; 
            t_4next = t;
            
            %Compute mutual invasibility
            %Set up resident as a mutant
            abun_vec_tmp =y(end,:)'; %Collect equilibrium population sizes 
            N_tot_tmp = reshape(abun_vec_tmp,m_tot,n); %Reshape abun_vec to abundanse matrix where rows are species and collumns are habitats 
            P_tmp=N_tot_tmp(m_cons+1:end,:);
            P_tmp(end+1,:)=[0,0,0]; %Ad a mutant (here the resident) 
            P_tmp(end,mutHab)=1; %Introduce one mutant individual in the correct habitat
            P_tmp(mut_pop,mutHab) = P_tmp(mut_pop,mutHab)-1; 
            Z_tmp(end+1)=Z(mut_pop);
            
            N_tmp = N_tot_tmp(1:m_cons,:);
            
            [m_pred_tmp n]=size(P_tmp); %no. rows represent no of species and no collumns represent no habitats
            N_tot_tmp=[N_tmp;P_tmp];
            abun_vec_tmp = reshape(N_tot_tmp,1,(m_cons+m_pred_tmp)*n)';

            parameters.N_tot = N_tot_tmp;
            parameters.N = N_tmp;
            parameters.P = P_tmp;
            parameters.abun_vec=abun_vec_tmp; %Collect equilibrium population sizes 
            parameters.mutHab=mutHab;
            parameters.Z=Z_tmp;

            %Compute resident invasion fitness            
            [t,y,mutFit_loc_pred mutFit_glob_pred mutFit_sim_pred mutInvFit_pred]=pred_fitness_func(parameters,100);
                        
             if mutFit_sim_pred > 0 %If mutual invasibility exist
               
                disp('Mutual invasibility exist') 
                %Compute new equilibrium with both resident and mutant present
                parameters.N_tot = [N; P_mutinv];
                parameters.N = N;
                parameters.P = P_mutinv;
                parameters.abun_vec=abun_vec_mutinv; %Collect equilibrium population sizes 
                parameters.Z=Z_mutinv;
                
                [t,y,K_a,alpha,a]=pop_dyn_func(parameters,10000);

                %Update data for next evolutionary step
                 %Save data for next evolutionary step
                Z_4next = Z_mutinv;
                pred_sp_id4next=[pred_sp_id pred_sp_id(mut_pop)];
                y(y<0) = 0;%Remove negative values, even though they should not occur according to the ODE solver settings
                y_4next = y; 
                t_4next = t;
             end
             
             %Update abundance and trait data
            Z = Z_4next;
            pred_sp_id=pred_sp_id4next;
            m_pred = length(Z);
            abun_vec=y_4next(end,:)'; %Collect equilibrium population sizes 
            m_tot=m_cons+m_pred;
            N_tot = reshape(abun_vec,m_tot,n); %Reshape abun_vec to abundanse matrix where rows are species and collumns are habitats 
            N=N_tot(1:m_cons,:);
            P=N_tot(m_cons+1:end,:);
            
            
            %Perge the system from populations with abudnance < 1
            [N V P Z m_cons n m_pred N_tot m_tot abun_vec prey_sp_id pred_sp_id]=PergeFunc(N,P,V,Z,prey_sp_id,pred_sp_id);
            
            
            %Recompute the rate w 
            tmp1 = N.*r*mut_N;
            tmp2 = P.*abs(d)*mut_P;
            tmp3 = [tmp1;tmp2];
            w = reshape(tmp3,1,m_tot*n)';
            w_tot=sum(w);

                        
        end % end of if invasion fitness is positive test
        
    end %end of if/ elseif loop
    
    %Check for speciation
    if check_PredPhylo_flag==1 | check_PreyPhylo_flag==1 
        [prey_sp_id, pred_sp_id, prey_sp_id_counter, pred_sp_id_counter, prey_phylo_data, pred_phylo_data]=PhyloCheck(V, prey_sp_id, Z, pred_sp_id,prey_sp_id_counter,pred_sp_id_counter, sigma_mut_N, sigma_mut_P, check_PreyPhylo_flag,check_PredPhylo_flag,t_evo, prey_phylo_data, pred_phylo_data);
    end    
    
    %Save data     
    prey_dist_data{end+1,1}=t_evo; %time
    prey_dist_data{end,2}=[V; prey_sp_id]; %matrix for trait distribution and species id
    prey_dist_data{end,3}=N; %Abundance for each morph in each habitat

    pred_dist_data{end,1}=t_evo; %time
    pred_dist_data{end,2}=[Z; pred_sp_id]; %matrix for trait distribution and species id
    pred_dist_data{end,3}=P; %Abundance for each morph in each habitat

    com_info_data{end+1,1}=t_evo; %time
    com_info_data{end,2}=K_a; %K_A
    com_info_data{end,3}=alpha; %alpha
    com_info_data{end,4}=a; %a
    
    
    %Do some plotting ones in a while
    if rem(t_evo,1)==0        
        PlotTraitPhylo(N,P,V,Z,t_evo)
    end %end of plot statement

end %end of evolutonary time while loop

save(name,'prey_dist_data','prey_phylo_data','pred_dist_data','pred_phylo_data','com_info_data');

