function [t y mutFit_loc mutFit_glob mutFit_sim mutInvFit] = fitness_func(parameters,t_end)



%This function computes fitness for a given mutant that appears in a given habitat

%Extract the parameters for the resident species
N = parameters.N;
P = parameters.P;
N_tot = parameters.N_tot;
abun_vec=parameters.abun_vec;
V=parameters.V; %Initial conditions traits
Z=parameters.Z; %Initial conditions traits
r=parameters.r; %Consumer growth
d=parameters.d; %Predator death
K0=parameters.K0; %Max K
bmax=parameters.bmax; %Max predator consumption
sigma_a=parameters.sigma_a; %Width consumption kernel
sigma_b=parameters.sigma_b; %Width predation kernel
mN=parameters.mN; %Consumer dispersal
mP=parameters.mP; %Predator dispersal
cP=parameters.cP; %Predator conversion rate
U=parameters.U; %Resource peaks
sigma_K=parameters.sigma_K; %Resource width
mutHab=parameters.mutHab; %Habitat where mutant appears

m_cons=length(V); %No consumers
m_pred = length(Z); %No predators
m_tot = m_cons+m_pred; %No species
n = length(K0); %No habitats
Vmut=V(end); %Mutant trait

%Set up the growht vector for mutant
growth_vec=ones(n,1)*r;

%Compute Carrying capacity for mutant in each habitat
K_vec=[];
for i=1:n %loop over habitats
        K_vec(i,1) = K0(i)*exp(-(Vmut-U(i)).^2/2/sigma_K(i)^2); %K for both species in habitat A
end

%%Set up the interaction between mutant and all resident consumers
alpha=[];
for i=1:m_cons
        alpha(1,i)=exp(-(V(i)-Vmut).^2/2/sigma_a^2);  
end

%%Set up the "from prey to predator" interaction
a=[]; 
for j=1:m_pred
    a(1,j)=-bmax*exp(-(Vmut-Z(j)).^2/2/sigma_b^2); 
end

%Put all interactions together in one matrix
%Note, here we multiply -r with alpha and we devide alpha with K to get the
%correct r-alpha-K formulation (see Case pp 347)
A=zeros(n,m_cons*n);
pos =1; 
for i =1:n
    A(i,pos:pos+m_tot-1)=[-r*alpha/K_vec(i) a];
    pos = pos + m_tot;
end

%Set up the dispersal matrix
disp_mat1=diag(ones(n,1)*-mN);
disp_mat2=diag(ones(n-1,1)*(mN),1);
disp_mat3=diag(ones(n-1,1)*(mN),-1);

disp_mat=disp_mat1+disp_mat2+disp_mat3;


%Mutant abundance vector
abun_vec_mut = N(end,:)';

y0=abun_vec_mut;

%Save the newly formed matrices
parameters.A = A; 
parameters.disp_mat = disp_mat; 
parameters.abun_vec = abun_vec; 
parameters.growth_vec = growth_vec;
parameters.abun_vec_mut = abun_vec_mut;

%Compute local fitness
tmp = abun_vec_mut.*(growth_vec+(A*abun_vec));
mutFit_loc = tmp(mutHab);

%Compute global fitness
tmp = abun_vec_mut.*(growth_vec+(A*abun_vec))+(disp_mat*abun_vec_mut);
mutFit_glob = sum(tmp);


%Compute simulated fitness
y0=abun_vec_mut; 
ode_opts = odeset('NonNegative',1:n); %options for the ode, do not let the population densities below zero

%Call the ode15 solver and input ode_sys function which contains the diff.
%equations

[t,y] = ode15s( @(t,y) ode_sys_mutfit_consumer( t , y , parameters ) , [0,t_end] , y0 , ode_opts ); %Here we use ode15s as it works better than ode45 with the non-negative options

% figure(222)
% plot( t , y)
% title('mutant initial growth')

if length(find(sum(y,2)<1))>0
    mutFit_sim = 0;
    mutInvFit=0;
else
    val = t_end-(t_end*0.3); %value to find
    tmp = abs(t-val);
    [idx idx] = min(tmp); %index of closest value
    closest = t(idx); %closest value

    mutFit_sim =(log(sum(y(end,:))/sum(y(idx,:))))/ (t(end)-closest);
    
    if mutFit_sim < 0.000001
        mutInvFit=0;
    elseif mutFit_sim > 0 
        mutInvFit = 1; 
    end
end