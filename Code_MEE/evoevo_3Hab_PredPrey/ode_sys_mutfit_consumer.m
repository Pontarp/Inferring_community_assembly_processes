function fy = ode_sys( t , y , parameters )

r=parameters.growth_vec;
N_mut = y;
N_res = parameters.abun_vec;
A = parameters.A; 
m = parameters.disp_mat;

dN=N_mut.*(r+(A*N_res))+(m*N_mut);

fy = [dN];