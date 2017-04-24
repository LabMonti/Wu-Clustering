function [index]=PrintCDbwPartition(X1,order,partitions)
prompt = {'Enter Number of Partition you Wish to Extract: ','Enter Title: ', 'Enter x-axis label: ','Enter y-axis label: '};
dlg_title = 'Input';
defaultans = {'1','Title','Displacement (nm)','Conductance (G/Go)'};
usrinput = inputdlg(prompt,dlg_title,1,defaultans);
index = str2num(usrinput{1});
usrtitle = usrinput{2};
usrx = usrinput{3};
usry = usrinput{4};

Y = partitions{index};

map = distinguishable_colors(max(Y(:))+2);
map1 = distinguishable_colors(max(Y(:)));
map(1,:) = [1 1 1];

figure();
scatter3(X1(order,1),10.^(X1(order,2)),X1(order,3), 50, Y, 'filled')
view(2)
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