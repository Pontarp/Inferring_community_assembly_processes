clear all
close all


load Extract_sig_a1rep1_model_2.mat

%Remove runs that did not qualify as a comparison agains the true data
%This can e.g. be if the predators crashed or if there is no prey diversity
tmp = find(isnan(SS_sim_mat(:,end))==1);
SS_sim_mat(tmp,:)=[];

%Convert and plot scatter data in surface
x=SS_sim_mat(:,1);
y=SS_sim_mat(:,2);
z=SS_sim_mat(:,end);

% Use Delaunay triangulation.
tri = delaunay(x,y);
[r,c] = size(tri); % How many triangles are there?
% Plot it with TRISURF
figure(1)
% subplot(3,1,1)
h = trisurf(tri, x, y, z);

% Clean up the plot
axis vis3d
l = light('Position',[-50 -15 29])
shading interp
colorbar EastOutside
xlabel('sig-alpha'); ylabel('sig_a'); zlabel('dist')
title('True = Pred (0.1) : Model = Pred-Prey')

%-----------------------------------
%Plot the euclidean distance for different sig_alpha, given that sig_a is correct
tmp = find(SS_sim_mat(:,2)>0.68 & SS_sim_mat(:,2)<0.72);
figure(11)
subplot(2,1,1)
plot(SS_sim_mat(tmp,1),SS_sim_mat(tmp,end),'r*');
axis square
ylabel('Euclidean distance')
xlabel('Sigma alpha')
title('True = Pred (0.1) : Model = Pred-Prey')

%Plot the euclidean distance for different sig_a, given that sig_alpha is correct
tmp = find(SS_sim_mat(:,1)>0.08 & SS_sim_mat(:,1)<0.11);
subplot(2,1,2)
plot(SS_sim_mat(tmp,2),SS_sim_mat(tmp,end),'r*');
axis square
ylabel('Euclidean distance')
xlabel('Sigma a')

%--------------------------------------------------
%Plott the data with marginal probability distribution
epsi=5%1.5*min(SS_sim_mat(:,end));

figure(111)
tmp=find(SS_sim_mat(:,end)<epsi);
y=SS_sim_mat(tmp,1);
x=SS_sim_mat(tmp,2);
scatterhist(x,y,'kernel','on','location','northeast','direction','out')
ylabel('sig alpha'); xlabel('sig a')


%%
clear all

%Load true data
load Extract_sig_a2rep1_model_2.mat

%Remove runs that did not qualify as a comparison agains the true data
%This can e.g. be if the predators crashed or if there is no prey diversity
tmp = find(isnan(SS_sim_mat(:,end))==1);
SS_sim_mat(tmp,:)=[];

%Convert and plot scatter data in surface
x=SS_sim_mat(:,1);
y=SS_sim_mat(:,2);
z=SS_sim_mat(:,end);

% Use Delaunay triangulation.
tri = delaunay(x,y);
[r,c] = size(tri); % How many triangles are there?
% Plot it with TRISURF
figure(2)
% subplot(3,1,1)
h = trisurf(tri, x, y, z);

% Clean up the plot
axis vis3d
l = light('Position',[-50 -15 29])
shading interp
colorbar EastOutside
xlabel('sig-alpha'); ylabel('sig_a'); zlabel('dist')
title('True = Pred (0.3) : Model = Pred-Prey')

%-----------------------------------
%Plot the euclidean distance for different sig_alpha, given that sig_a is correct
tmp = find(SS_sim_mat(:,2)>0.68 & SS_sim_mat(:,2)<0.72);
figure(22)
subplot(2,1,1)
plot(SS_sim_mat(tmp,1),SS_sim_mat(tmp,end),'r*');
axis square
ylabel('Euclidean distance')
xlabel('Sigma alpha')
title('True = Pred (0.3) : Model = Pred-Prey')


%Plot the euclidean distance for different sig_a, given that sig_alpha is correct
tmp = find(SS_sim_mat(:,1)>0.08 & SS_sim_mat(:,1)<0.11);
subplot(2,1,2)
plot(SS_sim_mat(tmp,2),SS_sim_mat(tmp,end),'r*');
axis square
ylabel('Euclidean distance')
xlabel('Sigma a')

%--------------------------------------------------
%Plott the data with marginal probability distribution
epsi=5%1.5*min(SS_sim_mat(:,end));

figure(222)
tmp=find(SS_sim_mat(:,end)<epsi);
y=SS_sim_mat(tmp,1);
x=SS_sim_mat(tmp,2);
scatterhist(x,y,'kernel','on','location','northeast','direction','out')
ylabel('sig alpha'); xlabel('sig a')




%%
clear all

load Extract_sig_a1sig_b4bmax3mutP2rep1_model_2.mat

%Remove runs that did not qualify as a comparison agains the true data
%This can e.g. be if the predators crashed or if there is no prey diversity
tmp = find(isnan(SS_sim_mat(:,end))==1);
SS_sim_mat(tmp,:)=[];

%Convert and plot scatter data in surface
x=SS_sim_mat(:,1);
y=SS_sim_mat(:,2);
z=SS_sim_mat(:,end);

% Use Delaunay triangulation.
tri = delaunay(x,y);
[r,c] = size(tri); % How many triangles are there?
% Plot it with TRISURF
figure(3)
% subplot(3,1,1)
h = trisurf(tri, x, y, z);

% Clean up the plot
axis vis3d
l = light('Position',[-50 -15 29])
shading interp
colorbar EastOutside
xlabel('sig-alpha'); ylabel('sig_a'); zlabel('dist')
title('True = Pred-Prey (0.7:0.1) : Model = Pred-Prey')

%-----------------------------------
%Plot the euclidean distance for different sig_alpha, given that sig_a is correct
tmp = find(SS_sim_mat(:,2)>0.68 & SS_sim_mat(:,2)<0.72);
figure(33)
subplot(2,1,1)
plot(SS_sim_mat(tmp,1),SS_sim_mat(tmp,end),'*');
axis square
ylabel('Euclidean distance')
xlabel('Sigma alpha')
title('True = Pred-Prey (0.7:0.1) : Model = Pred-Prey')

%Plot the euclidean distance for different sig_a, given that sig_alpha is correct
tmp = find(SS_sim_mat(:,1)>0.08 & SS_sim_mat(:,1)<0.11);
subplot(2,1,2)
plot(SS_sim_mat(tmp,2),SS_sim_mat(tmp,end),'*');
axis square
ylabel('Euclidean distance')
xlabel('Sigma a')

%----------------------------------------
%Plott the data with marginal probability distribution
epsi=1%1.5*min(SS_sim_mat(:,end));

figure(333)
tmp=find(SS_sim_mat(:,end)<epsi);
y=SS_sim_mat(tmp,1);
x=SS_sim_mat(tmp,2);
scatterhist(x,y,'kernel','on','location','northeast','direction','out')
ylabel('sig alpha'); xlabel('sig a')

%%
clear all

%Load true data
load Extract_sig_a2sig_b4bmax3mutP2rep1_model_2.mat
%Remove runs that did not qualify as a comparison agains the true data
%This can e.g. be if the predators crashed or if there is no prey diversity
tmp = find(isnan(SS_sim_mat(:,end))==1);
SS_sim_mat(tmp,:)=[];

%Convert and plot scatter data in surface
x=SS_sim_mat(:,1);
y=SS_sim_mat(:,2);
z=SS_sim_mat(:,end);

% Use Delaunay triangulation.
tri = delaunay(x,y);
[r,c] = size(tri); % How many triangles are there?
% Plot it with TRISURF
figure(4)
%subplot(3,1,1)
h = trisurf(tri, x, y, z);

% Clean up the plot
axis vis3d
l = light('Position',[-50 -15 29])
shading interp
colorbar EastOutside
xlabel('sig-alpha'); ylabel('sig_a'); zlabel('dist')
title('True = Pred-Prey (0.7:0.3) : Model = Pred-Prey')

%-----------------------------------
%Plot the euclidean distance for different sig_alpha, given that sig_a is correct
tmp = find(SS_sim_mat(:,2)>0.68 & SS_sim_mat(:,2)<0.72);
figure(44)
subplot(2,1,1)
plot(SS_sim_mat(tmp,1),SS_sim_mat(tmp,end),'*');
axis square
ylabel('Euclidean distance')
xlabel('Sigma alpha')
title('True = Pred-Prey (0.7:0.3) : Model = Pred-Prey')

%Plot the euclidean distance for different sig_a, given that sig_alpha is correct
tmp = find(SS_sim_mat(:,1)>0.08 & SS_sim_mat(:,1)<0.11);
subplot(2,1,2)
plot(SS_sim_mat(tmp,2),SS_sim_mat(tmp,end),'*');
axis square
ylabel('Euclidean distance')
xlabel('Sigma a')

%--------------------------------------------------
%Plott the data with marginal probability distribution
epsi=1%1.5*min(SS_sim_mat(:,end));

figure(444)
tmp=find(SS_sim_mat(:,end)<epsi);
y=SS_sim_mat(tmp,1);
x=SS_sim_mat(tmp,2);
scatterhist(x,y,'kernel','on','location','northeast','direction','out')
ylabel('sig alpha'); xlabel('sig a')