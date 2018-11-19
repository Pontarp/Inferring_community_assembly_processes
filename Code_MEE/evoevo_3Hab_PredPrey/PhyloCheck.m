function [prey_sp_id, pred_sp_id, prey_sp_id_counter, pred_sp_id_counter, prey_phylo_data, pred_phylo_data]=PhyloCheck(V, prey_sp_id, Z, pred_sp_id,prey_sp_id_counter,pred_sp_id_counter, sigma_mut_N, sigma_mut_P, check_PreyPhylo_flag,check_PredPhylo_flag,t_evo, prey_phylo_data, pred_phylo_data)

origin=[];
%Update the prey species id and phylogenetic inforamation

if check_PreyPhylo_flag ==1
    V_id_combo=[V; prey_sp_id]; %put trait distribution and species id together
    V_id_combo_sort=sortrows(V_id_combo')'; %Sort the matrix according to the trait distribuiton
    [dummy,p]=sort(V); %p provides the element indeces of V before sorted, will help us turn  the data back to original form

    sp_pres=unique(V_id_combo(2,:)); %Get present species to loop over

    for j=sp_pres %loop over species to check for gaps in them
        sp_j_ind=find(V_id_combo_sort(2,:)==j); %get index of specis j

        V_diff=abs(diff(V_id_combo_sort(1,sp_j_ind))); %compute the distance between each trait value in the j'th species trait distribution
        gap=find(V_diff > sigma_mut_N*3); %find gaps in V that are > x times the size of the mutations

        if length(gap)>0 %if gaps exist in the trait distribution vector

            tmp1=sp_j_ind(1:gap); %Get elements up to the gap, for the given species j
            tmp2=sp_j_ind(gap+1:end); %Get elements past the gap, for the given species j

            V_id_combo_sort(2,tmp1)=prey_sp_id_counter+1; %give the elements up to the gap a new id
            V_id_combo_sort(2,tmp2)=prey_sp_id_counter+2; %give the elements beyond the gap a new id

            prey_sp_id_counter=prey_sp_id_counter+2; %update the counter
            
            origin=j; 

            %Fill in phylo data 
            prey_phylo_data(end+1,1)=prey_sp_id_counter-1; %species id
            prey_phylo_data(end,2)=j; %origin
            prey_phylo_data(end,3)=t_evo; %time first registered  
            
            prey_phylo_data(end+1,1)=prey_sp_id_counter; %species id
            prey_phylo_data(end,2)=j; %origin
            prey_phylo_data(end,3)=t_evo; %time first registered  
            
        end
    end

    %Change the data back to original form
    tmp=[p; V_id_combo_sort];
    tmp2=sortrows(tmp',1)';
    prey_sp_id=tmp2(3,:); %prey_sp_id is the only vector that needs updating, N0 and Z should be the same

end

%%
%Update the pred species id and phylogenetic inforamation

if check_PredPhylo_flag==1
    Z_id_combo=[Z; pred_sp_id]; %put trait distribution and species id together
    Z_id_combo_sort=sortrows(Z_id_combo')'; %Sort the matrix according to the trait distribuiton
    [dummy,p]=sort(Z); %p provides the element indeces of V before sorted, will help us turn  the data back to original form

    sp_pres=unique(Z_id_combo(2,:)); %Get present species to loop over

    for j=sp_pres %loop over species to check for gaps in them
        sp_j_ind=find(Z_id_combo_sort(2,:)==j); %get index of specis j

        Z_diff=abs(diff(Z_id_combo_sort(1,sp_j_ind))); %compute the distance between each trait value in the j'th species trait distribution
        gap=find(Z_diff > sigma_mut_P*3); %find gaps in Z that are > x times the size of the mutations

        if length(gap)>0 %if gaps exist in the trait distribution vector

            tmp1=sp_j_ind(1:gap); %Get elements up to the gap, for the given species j
            tmp2=sp_j_ind(gap+1:end); %Get elements past the gap, for the given species j

            Z_id_combo_sort(2,tmp1)=pred_sp_id_counter+1; %give the elements up to the gap a new id
            Z_id_combo_sort(2,tmp2)=pred_sp_id_counter+2; %give the elements beyond the gap a new id

            pred_sp_id_counter=pred_sp_id_counter+2; %update the counter
            
            origin=j; 

            %Fill in phylo data 
            pred_phylo_data(end+1,1)=pred_sp_id_counter-1; %species id
            pred_phylo_data(end,2)=j; %origin
            pred_phylo_data(end,3)=t_evo; %time first registered
            
            pred_phylo_data(end+1,1)=pred_sp_id_counter; %species id
            pred_phylo_data(end,2)=j; %origin
            pred_phylo_data(end,3)=t_evo; %time first registered
        end
    end

    %Change the data back to original form
    tmp=[p; Z_id_combo_sort];
    tmp2=sortrows(tmp',1)';
    pred_sp_id=tmp2(3,:); %pred_sp_id is the only vector that needs updating, P0 and V should be the same
end
    
    
    
    