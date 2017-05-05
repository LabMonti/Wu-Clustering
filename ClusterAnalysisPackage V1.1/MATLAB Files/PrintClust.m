function PrintClust(X1,Y,usrtitle,usrx,usry)

map = distinguishable_colors(max(Y(:))+2);
map1 = distinguishable_colors(max(Y(:)));
map(1,:) = [1 1 1];
map(2,:) = [1 1 1];


figure();
scatter3(X1(:,1),10.^(X1(:,2)),X1(:,3), 50, Y, 'filled')
view(2)
set(gca,'yscale','log')
grid off

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
set(gca,'fontsize',16)
title(usrtitle,'fontweight','bold','fontsize',24)
xlabel(usrx,'fontsize',16) % x-axis label
ylabel(usry,'fontsize',16)
end

