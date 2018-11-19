function PlotInitFitLand(mut_range,mut_hab,abun_vec,N_tot,N,P,V,Z,parameters)

%Extract parameters
m_pred = parameters.m_pred;
m_cons = parameters.m_cons;
n=parameters.n; 


%Variables for fitness values
mutFit_loc_cons=[];
mutFit_glob_cons=[];
mutFit_sim_cons=[];
mutInvFit_cons = [];
mutFit_loc_pred=[];
mutFit_glob_pred=[];
mutFit_sim_pred=[];
mutInvFit_pred = [];

for i=mut_hab
    
    mutHab=i; %The habitat where the species occurs
    mut_pop = i; %The species that mutated, here we assume that the best adapted species is mutating
    
    %-------When consumer mutate---------------
    N_tmp = N; 
    N_tmp(end+1,:)=[0,0,0]; %Ad a mutant 
    N_tmp(end,mutHab)=1; %Introduce one mutant individual in the correct habitat
    N_tmp(mut_pop,mutHab) = N_tmp(mut_pop,mutHab)-1; %Remove a individual from the mutating species
    V_tmp = V; %Set up trait vector that includes mutant
    V_tmp(end+1)=0; 
   
    [m_cons_tmp n]=size(N_tmp); %no rows represent no of species and no collumns represent no habitats
    N_tot_tmp=[N_tmp;P];
    abun_vec = reshape(N_tot_tmp,1,(m_cons_tmp+m_pred)*n)';
    
    parameters.N_tot = N_tot_tmp
    parameters.N = N_tmp;
    parameters.P = P;
    parameters.Z=Z;
    parameters.abun_vec=abun_vec; %Collect equilibrium population sizes 
    parameters.mutHab=mutHab;

    count=0;
        for j = mut_range;
            count=count+1;
            
            V_tmp(end)=j;
            parameters.V=V_tmp;

            %Compute invasion fitness            
            [t,y,mutFit_loc_cons(i,count) mutFit_glob_cons(i,count) mutFit_sim_cons(i,count) mutInvFit_cons(i,count)]=cons_fitness_func(parameters,100);            
        end
        
        
      %-------When predators mutate---------------
    P_tmp = P; 
    P_tmp(end+1,:)=[0,0,0]; %Ad a mutant 
    P_tmp(end,mutHab)=1; %Introduce one mutant individual in the correct habitat
    P_tmp(mut_pop,mutHab) = P_tmp(mut_pop,mutHab)-1; %Remove a individual from the mutating species
    Z_tmp = Z; %Set up trait vector that includes mutant
    Z_tmp(end+1)=0;
 
    [m_pred_tmp n]=size(P_tmp); %no rows represent no of species and no collumns represent no habitats
    N_tot_tmp=[N;P_tmp];
    abun_vec = reshape(N_tot_tmp,1,(m_cons+m_pred_tmp)*n)';
    
    parameters.N_tot = N_tot_tmp;
    parameters.N = N;
    parameters.V = V; 
    parameters.P = P_tmp;
    parameters.abun_vec=abun_vec; %Collect equilibrium population sizes 
    parameters.mutHab=mutHab;

    count=0;
        for j = mut_range;
            count=count+1;
            
            Z_tmp(end)=j;
            parameters.Z=Z_tmp;
            
            %Compute invasion fitness
            [t,y,mutFit_loc_pred(i,count) mutFit_glob_pred(i,count) mutFit_sim_pred(i,count) mutInvFit_pred(i,count)]=pred_fitness_func(parameters,100);

        end  
end

%Plot fitness landscape
figure(4)
subplot(2,4,1)
plot(mut_range,mutFit_loc_cons'); ylim([-0.5 1]); legend('Hab1','Hab2','Hab3')
title('Fitness locally'); ylabel('Fitness'); xlabel('Trait')
subplot(2,4,2)
plot(mut_range,mutFit_glob_cons'); ylim([-0.5 1]);
title('Fitness Globally'); xlabel('Trait')
subplot(2,4,3)
plot(mut_range,mutFit_sim_cons');
title('Fitness Simulated'); xlabel('Trait')
subplot(2,4,4)
plot(mut_range, mutInvFit_cons'); ylim([-0.5 1.5]);
title('Invasion probability'); xlabel('Trait')
subplot(2,4,5)
plot(mut_range,mutFit_loc_pred'); ylim([-0.5 1]); legend('Hab1','Hab2','Hab3')
title('Fitness locally'); ylabel('Fitness'); xlabel('Trait')
subplot(2,4,6)
plot(mut_range,mutFit_glob_pred'); ylim([-0.5 1]);
title('Fitness Globally'); xlabel('Trait')
subplot(2,4,7)
plot(mut_range,mutFit_sim_pred');
title('Fitness Simulated'); xlabel('Trait')
subplot(2,4,8)
plot(mut_range, mutInvFit_pred'); ylim([-0.5 1.5]);
title('Invasion probability'); xlabel('Trait')
