function Print2DHist(X1,usrtitle,usrx,usry)
load('cmap.mat');
figure()
scatter(X1(:,1),10.^(X1(:,2)),50,X1(:,3),'fill','o')
set(gca,'yscale','log')
set(gca,'fontsize',16)
xlabel(usrx,'fontsize',16) % x-axis label
ylabel(usry,'fontsize',16) % y-axis label
title(usrtitle,'fontweight','bold','fontsize',24)
colormap(cmap)
colorbar
end