function [partition,index] = CDbwAnalysis(X1,X2,RD,CD,order,geneps,eps,n)

h = waitbar(0,'Please wait...');

%Extract a partition (clustering solution) for each value in the vector eps
%Store cluster ID vector and CDbw value for each partition calculated
for i=1:n
    [Y1,~,CDbw] = ExtractOpticsPartition(X1,X2,RD,CD,order,eps(i),geneps);
    drawnow
    partition{i} = Y1;
    index(i) = CDbw;
    
    waitbar(i/n);
end

close(h);

figure()
plot(eps,index,'-o','markerfacecolor','k')
set(gca,'fontsize',16)
xlabel('Epsilon')
ylabel('CDbw')
title('Determining Optimal Partition','fontweight','bold','fontsize',24)

[~,I] = max(index);

Y = partition{I};
figure()
scatter3(X1(order,1),10.^(X1(order,2)),X1(order,3), 50, Y, 'filled')
view(2)
set(gca,'yscale','log')
grid off

map = distinguishable_colors(max(Y(:))+2);
map1 = distinguishable_colors(max(Y(:)));
map(1,:) = [1 1 1];
map(2,:) = [1 1 1];

if min(Y(:)) == -1
    caxis(caxis)
    colormap(map) 
    labels(1) = -1;
    for l = 2:(max(Y(:))+2)
        labels(l) = l-2;
    end
else
    colormap(map1) 
    for l = 1:(max(Y(:)))
        labels(l) = l;
    end
end
lcolorbar(labels,'FontSize', 14,'fontweight','bold');
title('Optimal CDbw Partition','FontWeight','bold','fontsize',18)
ylabel('Conductance (G/Go)','Fontweight','bold','fontsize',14)
xlabel('Distance (m)','FontWeight','bold','fontsize',14) % x-axis label
set(gca,'fontweight','bold')


end