function TraitPhyloPlot(N,P,V,Z,t_evo)

%Do some plotting ones in a while

 %Plot Consumer trait distributions and abundance--------------
        %rgb color coding according to the abudnance in different habitats
        [mm nn]=size(N);
        c=[];
        for i=1:mm
            c(i,:)=N(i,:)/sum(N(i,:));
        end
        
        figure(200)
        subplot(2,2,1)
        cla
        hold on
        for i=1:mm
%            bar(V(i),sum(N(i,:)),'barwidth',0.4,'facecolor',c(i,:),'edgecolor','none','facealpha',0.7); %only works in latest version of Matlab
            bar(V(i),sum(N(i,:)),'barwidth',0.4,'facecolor',c(i,:),'edgecolor','none');%works also in MATLAB 2014 
        end
        axis square
        subplot(2,2,2)
        hold on
        scatter(V,ones(1,mm)*t_evo,[],c,'fill','o');
        axis square
        
        if sum(sum(P))>1

            %Plot Predator trait distributions and abundance--------------
            %rgb color coding according to the abudnance in different habitats
            [mm nn]=size(P);
            c=[];
            for i=1:mm
                c(i,:)=P(i,:)/sum(P(i,:));
            end

            %Plot trait distribution and radiation
            figure(200)
            subplot(2,2,3)
            cla
            hold on
            for i=1:mm
    %             bar(Z(i),sum(P(i,:)),'barwidth',0.4,'facecolor',c(i,:),'edgecolor','none','facealpha',0.7);
                bar(Z(i),sum(P(i,:)),'barwidth',0.4,'facecolor',c(i,:),'edgecolor','none'); %works also in MATLAB 2014
            end
            axis square
            subplot(2,2,4)
            hold on
            scatter(Z,ones(1,mm)*t_evo,[],c,'fill','o');
            axis square
        end %end of test for Predator abundance above 1