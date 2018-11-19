function PlotSpecDyn(y,t,abun_vec, parameters)

m_tot=parameters.m_tot;
m_cons = parameters.m_cons;
n = parameters.n; 

[mm nn] = size(y);
points=[];
figure(2)
subplot(2,1,2)
axis square
hold on
for i = 1:mm %loop over rows in y as rows represent time
    abun_vec=y(i,:)'; %Collect population sizes 
    N_tot = reshape(abun_vec,m_tot,n); %Reshape abun_vec to abundanse matrix where rows are species and collumns are habitats 
    N=N_tot(1:m_cons,:);
    P=N_tot(m_cons+1:end,:);
    
    %Plot prey abundances with rgb color coding according to the abudnance in different habitats
    tmp=sum(N,2);
    tmp2=find(tmp>0);
    for j=tmp2'
        c=N(j,:)/sum(N(j,:));
        plot(t(i),sum(N(j,:)),'-o','markerfacecolor',c,'markeredgecolor','none')
        points(j,i)=sum(N(j,:));
    end
    
    %Plot pred abundances with rgb color coding according to the abudnance in different habitats
    tmp=sum(P,2);
    tmp2=find(tmp>0);
    for j=tmp2'
        c=P(j,:)/sum(P(j,:));
        plot(t(i),sum(P(j,:)),'-<','markerfacecolor',c,'markeredgecolor','none')
        points(j,i)=sum(P(j,:));
    end
end
