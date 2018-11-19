function fy = ode_sys( t , y , parameters )

d=parameters.death_vec;
P_mut = y;
N_res = parameters.abun_vec_cons;
A = parameters.A; 
m = parameters.disp_mat;

dN=P_mut.*(d+(A*N_res))+(m*P_mut);

fy = [dN];