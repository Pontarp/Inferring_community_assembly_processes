function fy = ode_sys( t , y , parameters )

r=parameters.growdeath_vec;
N = y;
A = parameters.A; 
m = parameters.disp_mat; 

dN = N.*(r+(A*N))+(m*N);
% dN=N.*(r+((-r.*1./K).*(A*N)))+(m*N);

fy = [dN];