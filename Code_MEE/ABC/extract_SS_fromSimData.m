%This script extracts summary statistics for each simulation of a community
%with sigma_a and sigma_alpha pulled randomly from a normal distribution
%between 0-1. 

%The output will be a matrix containing the parameters (first two collumns)
%of the simulation, followed by values on: 

% No. prey
% Mean abundance of prey
% Width of the prey trait distribution
% Mean trait distance MTD
% mean nearest trait distance MNTD
% Mean phylogenetic distance
% Nearest neighbor phylogenetic distance

close all
clear all

%Specify evolutionary time which defines the end community
endtime=3000; 

%%
%First collect SS from "true" data
%;Varable for ture data output
SS_true=[];

 load true_mainout_sig_a1rep1.mat
name=['Extract_sig_a1rep1_model_2'];

%  load true_mainout_sig_a2rep1.mat
% name=['Extract_sig_a2rep1_model_2'];


% load true_mainout_sig_a1sig_b4bmax3mutP2rep1.mat
% name=['Extract_sig_a1sig_b4bmax3mutP2rep1_model_2'];

%  load true_mainout_sig_a2sig_b4bmax3mutP2rep1.mat
% name=['Extract_sig_a2sig_b4bmax3mutP2rep1_model_2'];



%Get no. prey species in true data
tmp = prey_dist_data{endtime,2};
SS_true(1)=length(unique(tmp(3,:)));
%Get mean prey abundance
SS_true(2)=mean(tmp(2,:));
%Get width of the trait dist
SS_true(3)=abs(max(tmp(1,:))-min(tmp(1,:)));

%Pipe the community and phylo-info into the function that gives
%distance matrix
end_community=tmp;
[pair_phydist_mat pair_phydist_vec pair_traitdist_mat pair_traitdist_vec] = phyloconstruct(end_community, prey_phylo_data, endtime, 0);

%Pipe the distance matrix and vector into the function that
%gives MPD, NNMPD, MTD and MNTD
[MPD NNPD MTD MNTD] = speciesdiff(pair_phydist_mat, pair_phydist_vec, pair_traitdist_mat, pair_traitdist_vec);

SS_true(4)=MTD; 
SS_true(5)=MNTD;
SS_true(6)=MPD;
SS_true(7)=NNPD;

clearvars -except SS_true endtime name


%%
%Now collect SS from simulated data
%Variable for simulated data output
SS_sim_mat = [];

%Load simulated data
list=ls('mainout*');
[m n]=size(list);

count=0;
for i=1:m %loop over output files
    load(list(i,:)); %load i'th file nam
    
    [mm nn]=size(pred_dist_data);
    if mm>=endtime
    
    %Check if predators exist
    %If no predators do not continue
    tmp = pred_dist_data{endtime,2};
    pred_div=length(unique(tmp(3,:)));
    pred_abun=mean(tmp(2,:));
   
    
    if pred_div>0 & pred_abun>0
        count=count+1;
       %Get sigma parameters
        SS_sim_mat(count,1)=sigma_a;
        SS_sim_mat(count,2)=sigma_b; %this is where the predator nich width will be


        %Get no. prey species
        tmp = prey_dist_data{endtime,2};
        SS_sim_mat(count,3)=length(unique(tmp(3,:)));

        %Get mean prey abundance
        SS_sim_mat(count,4)=mean(tmp(2,:));

        %Get width of the trait dist
        SS_sim_mat(count,5)=abs(max(tmp(1,:))-min(tmp(1,:)));

        if SS_sim_mat(count,3)>1 %only do phylo and trait evaluation if there is more than one species
            %Pipe the community and phylo-info into the function that gives
            %distance matrix
            end_community=tmp;
            [pair_phydist_mat pair_phydist_vec pair_traitdist_mat pair_traitdist_vec] = phyloconstruct(end_community, prey_phylo_data, endtime, 0);

            %Pipe the distance matrix and vector into the function that
            %gives MPD, NNMPD, MTD and MNTD
            [MPD NNPD MTD MNTD] = speciesdiff(pair_phydist_mat, pair_phydist_vec, pair_traitdist_mat, pair_traitdist_vec);
            SS_sim_mat(count,6)=MTD; 
            SS_sim_mat(count,7)=MNTD;
            SS_sim_mat(count,8)=MPD;
            SS_sim_mat(count,9)=NNPD;
        else 
           SS_sim_mat(count,6)=NaN;
           SS_sim_mat(count,7)=NaN;
           SS_sim_mat(count,8)=NaN;
           SS_sim_mat(count,9)=NaN;       
        end

        clearvars -except SS_true SS_sim_mat n m list endtime count name
    end %end of if statement checking for predators
end
end %end of output file loop

%% Now compute the distance between each simulated data set and the true data

%Get standard deviation for each SS in the simulated data. We use this to
%scale the distance between SS_true and SS_simulated. Otherwise the SS that
%take large values and large std will be owerweighted in the eucledian
%distance analysis

SS_std = nanstd(SS_sim_mat(:,3:9)); 

[m n]=size(SS_sim_mat);
dist=[];
for i = 1:m    
    %Compute eucledian distance        
    p_dist_SS= abs(SS_true-SS_sim_mat(i,3:9))./SS_std;
    dist(i)=sqrt(sum(p_dist_SS).^ 2); %Eucledian distance between summary stats, absolut value, used for, threashold evaluation
end

SS_sim_mat = [SS_sim_mat dist']; %ad distances as a final collumn in SS_sim_mat
SS_sim_mat = sortrows(SS_sim_mat,1); %Sort whole matrix according to sigma

%Remove runs that did not qualify as a comparison agains the true data
%This can e.g. be if the predators crashed or if there is no prey diversity
tmp = find(isnan(SS_sim_mat(:,end))==1);
SS_sim_mat(tmp,:)=[];

%Plot 3D scatterplot
figure(1)
scatter3(SS_sim_mat(:,1),SS_sim_mat(:,2),SS_sim_mat(:,end))  

%Convert and plot scatter data in surface
x=SS_sim_mat(:,1);
y=SS_sim_mat(:,2);
z=SS_sim_mat(:,end);


% Use Delaunay triangulation.
tri = delaunay(x,y);
[r,c] = size(tri); % How many triangles are there?
% Plot it with TRISURF
figure(2)
h = trisurf(tri, x, y, z);

% Clean up the plot
axis vis3d
l = light('Position',[-50 -15 29])
shading interp
colorbar EastOutside
xlabel('sig-alpha'); ylabel('sig_a'); zlabel('dist')

%Plott posterior distribution given som treshold epsi
epsi=7%1.5*min(SS_sim_mat(:,end));

tmp=find(SS_sim_mat(:,end)<epsi);

figure(3)
subplot(2,1,1)
histogram(SS_sim_mat(tmp,1));
title('Marginal dist sig-alpha')
subplot(2,1,2)
histogram(SS_sim_mat(tmp,2));
title('Marginal dist sig-a')


%Save output
save([name '.mat'],'SS_sim_mat');

%Save figures
for i=[1 2 3]
   saveas(i,[name '_fig_' num2str(i)])
end
    
    
    
    
    
    
    
    
    
    
    
    
    
    