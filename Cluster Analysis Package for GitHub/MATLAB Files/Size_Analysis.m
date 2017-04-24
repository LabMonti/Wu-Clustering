% Cluster Size analysis

% Intialize Cluster Size Array
function [T,ClusterSize] = Size_Analysis(Y)
ClusterSize = zeros(max(Y)+1,3);
ClusterSize(1,1) = -1;

for i = 2:max(Y)+1
    ClusterSize(i,1) = i-1;
end

% size of noise cluster
Unresolved = -1; 
CU = find(Y==Unresolved);
sizeCU = length(CU);
percentCU = sizeCU/length(Y);
ClusterSize(1,2) = sizeCU;
ClusterSize(1,3) = percentCU*100;

k = 2; % used to index through ClusterSize

for i=1:max(Y)
   C = find(Y==i);
   sizeC(i) = length(C);
   percentC(i) = sizeC(i)/length(Y)*100;
   ClusterSize(k,2) = sizeC(i);
   ClusterSize(k,3) = percentC(i);
   k = k+1;
end

[~, order] = sort(ClusterSize(:,3),'descend');
ClusterSize = ClusterSize(order,:);

T=array2table(ClusterSize, 'VariableNames',{'Index', 'Num_Points', 'Percent_Total'});
end




