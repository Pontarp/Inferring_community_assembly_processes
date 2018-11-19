function [N V P Z m_cons n m_pred N_tot m_tot abun_vec prey_sp_id pred_sp_id]=PergeFunc(N,P,V,Z,prey_sp_id,pred_sp_id)

%Perge the system from populations with abudnance < 1

ext_ind = find(sum(N,2)<1);
N(ext_ind,:)=[];
V(ext_ind)=[];
prey_sp_id(ext_ind)=[];
[m_cons n]=size(N);

ext_ind = find(sum(P,2)<1);
P(ext_ind,:)=[];
Z(ext_ind)=[];
pred_sp_id(ext_ind)=[];
[m_pred n]=size(P);

N_tot=[N;P];
[m_tot dumy]=size(N_tot);
abun_vec = reshape(N_tot,1,m_tot*n)';