function [RD,CD,order] = OpticsAnalysis(X,geneps,k)

[RD,CD,order] = opticsv2(X,geneps,k);

ax3 = subplot(2,2,1);
bar(sort(order),RD,1);

title('Reachability Plot','FontWeight','bold','fontsize',18)
xlabel('Clustering Order','FontWeight','bold','fontsize',14) % x-axis label
ylabel('Reachability Distance','FontWeight','bold','fontsize',14) % y-axis label
set(gca,'fontweight','bold')
colormap(ax3,jet);

end