function [pair_phydist_mat pair_phydist_vec pair_traitdist_mat pair_traitdist_vec] = phyloconstruct(end_community, phyloinfo, fulltime, treeflag)


%Define the end community, basically set which species to
%construct the tree for
sp_id=unique(end_community(3,:));

%Get evolutionary history for each species
H={}; %Cell for history, a vector with history for each species in each element of H
c1=0; %counter
for i=sp_id; %loop over species
    c1=c1+1; 
    h=NaN; %history vector for the focal species
    sp=i;
    c=0; %counter
    while h(end)~=0
        c=c+1;
        tmp=find(phyloinfo(:,1)==sp);
        
        h(c)=phyloinfo(tmp,2);
        sp=h(c);
    end    
    H{end+1}=h;
    
    names{c1,1}=['Species ' num2str(i)];
end

%Compute pairwise phylogenetic and phenotypic distance between species
D=[];
m=length(sp_id);
for i=1:m %loop over species
    for j=1:m %loop over species
        if i~=j
            
            %Phylogenetic distance
            common=max(intersect(cell2mat(H(i)),cell2mat(H(j)))); %Find closest common ancestor
            
            tmp=find(phyloinfo(:,2)==common); %get index in phylo data for common ancestor
            
            time=phyloinfo(tmp(end),3); %Get the evolutionary time for the common split
            
            D(i,j)=fulltime-time; %Save distance in matrix
            
            %Phenotypic distance
            tmp=find(end_community(3,:)==sp_id(i)); %index of species i
            trait1=mean(end_community(1,tmp)); %define the trait of the species as the mean of all morphs
            tmp=find(end_community(3,:)==sp_id(j)); %index of species j
            trait2=mean(end_community(1,tmp)); %define the trait of the species as the mean of all morphs
            
            D_trait(i,j)=abs(trait1-trait2); %Save distance in matrix
            
        end
    end
end

%Set up distance for output
pair_phydist_mat=D; 
pair_traitdist_mat=D_trait; 

%Convert matrix to vector, this will be the input for the tree drawing
%funciton
[m n]=size(pair_phydist_mat);

pair_phydist_vec=[];
pair_traitdist_vec=[];
for i=1:n
    tmp=find(pair_phydist_mat(i,:)==0);
    pair_phydist_vec=[pair_phydist_vec pair_phydist_mat(i,tmp+1:end)];
    
    tmp=find(pair_traitdist_mat(i,:)==0);
    pair_traitdist_vec=[pair_traitdist_vec pair_traitdist_mat(i,tmp+1:end)];
end

% names{1,1}='test1';
% names{2,1}='test2';
% names{3,1}='test3';
% names{4,1}='test4';

if treeflag==1; 
    %Plot the phylogenetic tree
     tree=seqlinkage(pair_phydist_vec,'average',names)
     plot(tree,'orient','bottom')
end

            
            