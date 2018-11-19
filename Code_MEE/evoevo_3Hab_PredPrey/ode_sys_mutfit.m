function fy = ode_sys( t , y , parameters )

r=parameters.growth_vec;
N_mut = y;
N_res = parameters.abun_vec;
K = parameters.K_vec; 
A = parameters.alpha_mat; 
m = parameters.disp_mat;

dN=N_mut.*(r+((-r.*1./K).*(A*N_res)))+(m*N_mut);

fy = [dN];