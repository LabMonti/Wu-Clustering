function ExtractClust(X,Y,order,RD)

prompt = {'Enter Cluster Indices You Would Like to Extract As A Comma Separated List:','Enter Title: ', 'Enter x-axis label: ','Enter y-axis label: '};
dlg_title = 'Input';
defaultans = {'1','Extracted Clustering','Displacement (nm)','Conductance (G/Go)'};
usrinput = inputdlg(prompt,dlg_title,1,defaultans);
idx = str2num(usrinput{1});
usrtitle = usrinput{2};
usrx = usrinput{3};
usry = usrinput{4};

X = X(order,:);

X1 = [];
Y1 = [];
RD1 = [];
order1 = [];

for i = 1:length(idx)
    X_idx = X(Y == idx(i),:);
    RD_idx = RD(Y==idx(i),:);
    X1 = [X1 ; X_idx];
    Y1 = [Y1; idx(i).*ones(length(X_idx),1)];
    RD1 = [RD1; RD_idx];
    order1 = [order1;find(Y==idx(i))];
    labels(i) = idx(i);
end

if min(Y(:)) == -1
    map = distinguishable_colors(max(Y(:))+2);
    fig = figure();
    colormap(fig,map(labels+2,:))
    scatter3(X1(:,1),10.^(X1(:,2)),X1(:,3), 50, Y1, 'filled')
    view(2)
    set(gca,'yscale','log')
    lcolorbar(labels,'FontSize', 14,'fontweight','bold');
    set(gca,'fontsize',16)
    title(usrtitle,'fontweight','bold','fontsize',24)
    xlabel(usrx,'fontsize',16) % x-axis label
    ylabel(usry,'fontsize',16)
    grid off
else
    map = distinguishable_colors(max(Y(:)));
    fig = figure();
    colormap(fig,map(labels,:))
    scatter3(X1(:,1),10.^(X1(:,2)),X1(:,3), 50, Y1, 'filled')
    view(2)
    set(gca,'yscale','log')
    lcolorbar(labels,'FontSize', 14,'fontweight','bold');
    grid off
    set(gca,'fontsize',16)
    title(usrtitle,'fontweight','bold','fontsize',24)
    xlabel(usrx,'fontsize',16) % x-axis label
    ylabel(usry,'fontsize',16)
end

figure();
bar(order1,RD1,1);
set(gca,'fontsize',16)
title('Extracted Reachability Plot','FontWeight','bold','fontsize',24)
xlabel('Clustering Order','fontsize',16) % x-axis label
ylabel('Reachability Distance','fontsize',16) % y-axis label
colormap(jet);
    

end