function PrintRP(order,RD,usrtitle)
figure();
bar(sort(order),RD,1);
set(gca,'fontsize',16)
title(usrtitle,'FontWeight','bold','fontsize',24)
xlabel('Clustering Order','fontsize',16) % x-axis label
ylabel('Reachability Distance','fontsize',16) % y-axis label
colormap(jet);
end

