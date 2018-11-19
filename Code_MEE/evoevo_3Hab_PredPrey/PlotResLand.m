function PlotResLand(parameters)

%Extract parameters
V0=parameters.V;
Z0=parameters.Z;
K0=parameters.K0;
sigma_a=parameters.sigma_a;
sigma_b=parameters.sigma_b;
U=parameters.U;
sigma_K=parameters.sigma_K;
n=parameters.n;

%Plot the resource landscape and species niches for the initial conditions
count=0; 
range=-3:0.01:6;
for i=range
    count=count+1;
    for j=1:n
        K_A(count,j) = K0(j)*exp(-(i-U(j)).^2/2/sigma_K(j)^2);
    end
end

c1=0;
for i = range
    c1=c1+1;    
    c2=0; 
    for j = V0
        c2=c2+1;        
        N_niche(c1,c2)= 5000*exp(-(i-j).^2/2/sigma_a^2);        
    end 
end

c1=0;
for i = range
    c1=c1+1;    
    c2=0; 
    for j = Z0
        c2=c2+1;        
        P_niche(c1,c2)= 5000*exp(-(i-j).^2/2/sigma_b^2);        
    end 
end

figure(1)
plot_color='rgb'
hold on
for i = 1:length(K0)
    plot(range,K_A(:,i),plot_color(i));
end
for i=1:length(V0)
    plot(range,N_niche(:,i),['--' plot_color(i)])
end
for i=1:length(Z0)
    plot(range,P_niche(:,i),[':' plot_color(i)])
end
title('resource distribution and species niche')
axis square