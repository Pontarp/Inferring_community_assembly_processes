function [t y K_A alpha a] = pop_dyn_func(parameters,t_end)

%This function computes population dynamics for a given time frame

%Extract the parameters that were piped into this script
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

m_cons=length(V); %No consumers
m_pred=length(Z); %No predators
m_tot = m_cons+m_pred; %No total species
n = length(K0); %No habitats

%Set up the growht/death vector
growth_mat=ones(m_cons,n)*r;
death_mat=ones(m_pred,n)*d;
tmp=[growth_mat; death_mat];
growdeath_vec = reshape(tmp,1,m_tot*n)';

%Compute Carrying capacity for each consumer population in each habitat
K_A = [];
for i=1:n %loop over habitats
    for j=1:m_cons %loop over 
        K_A(j,i) = K0(i)*exp(-(V(j)-U(i)).^2/2/sigma_K(i)^2); %K for both species in habitat A
    end    
end

%%Set up the consumer interaction matrix
alpha=[];
for i=1:m_cons
    for j=1:m_cons
        alpha(i,j)=exp(-(V(i)-V(j)).^2/2/sigma_a^2);  
    end
end


%%Set up the "from prey to predator" matrix
a=[]; 
for i=1:m_cons
    for j=1:m_pred
        a(i,j)=-bmax*exp(-(V(i)-Z(j)).^2/2/sigma_b^2); 
    end
end

%%Set up the to predator from prey matrix
a_tmp=a.*-1; %change signe on a
a_tmp=a_tmp'; %use the transpose of a
ka=cP*a_tmp;

%Set up the matrices of the GLV (see pp 348 in Case)
%See also Case pp 387 for an example of the full model with dispersal
%Note that each of the matrices will be devided in subblocks here as we
%have two habitats. 

%Set up the r-alpha-K subblocks for the A matrix, one per habitat
for i = 1:n %loop over habitats
    
    K_tmp=repmat(K_A(:,i),1,m_cons); %Carrying capacty for each species in a given habitat as a matrix
    r_alpha_K(:,:,i)=-r*(alpha./K_tmp); %r_alpha_K for a given habitat
    
end

%Set up the full A matrix
A=[];
for i=1:n
    A=blkdiag(A,[r_alpha_K(:,:,i) a; ka zeros(m_pred)]);
end

%Set up the dispersal matrix
[mm, nn]=size(A);

disp_mat1=diag(ones(mm,1)*-mN);
disp_mat2=diag(ones(mm-m_tot,1)*(mN),m_tot);
disp_mat3=diag(ones(mm-m_tot,1)*(mN),-m_tot);

disp_mat=disp_mat1+disp_mat2+disp_mat3;

y0=abun_vec;

%Save the newly formed matrices
parameters.A = A; 
parameters.disp_mat = disp_mat; 
parameters.abun_vec = abun_vec; 
parameters.growdeath_vec = growdeath_vec; 


%Solve for the population dyanmics 

ode_opts = odeset('NonNegative',1:length(mm)); %options for the ode, do not let the population densities below zero

%Call the ode15 solver and input ode_sys function which contains the diff.
%equations

[t,y] = ode15s( @(t,y) ode_sys_popeq( t , y , parameters ) , [0,t_end] , y0 , ode_opts ); %Here we use ode15s as it works better than ode45 with the non-negative options
