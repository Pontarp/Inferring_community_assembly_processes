function [MPD NNPD MTD MNTD] = speciesdiff(pair_phydist_mat, pair_phydist_vec, pair_traitdist_mat, pair_traitdist_vec)

%Compute MPD
MPD=mean(pair_phydist_vec); %mean of the distance between all species

%Compute NNPD
[mmm nnn]=size(pair_phydist_mat);
tmp=0;
for iii=1:mmm %loop over rows in the distance matrix
    dist_sort=sort(pair_phydist_mat(iii,:)); %Sort the distances
    tmp=tmp+dist_sort(2); %The first value will be zero, the second value will be the distance to the nearest neighbour
end
tmp=tmp/mmm; %NNPD for this run is the mean of all nearest neighbour distances 
NNPD=tmp;

%Compute MTD
MTD=mean(pair_traitdist_vec); %mean of the distance between all species

%Compute MNTD
[mmm nnn]=size(pair_traitdist_mat);
tmp=0;
for iii=1:mmm %loop over rows in the distance matrix
    dist_sort=sort(pair_traitdist_mat(iii,:)); %Sort the distances
    tmp=tmp+dist_sort(2); %The first value will be zero, the second value will be the distance to the nearest neighbour
end
tmp=tmp/mmm; %NNPD for this run is the mean of all nearest neighbour distances 
MNTD=tmp;