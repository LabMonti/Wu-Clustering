% Extracting DBSCAN clusters
% Algorithm based off of the paper "OPTICS: Ordering Points To Identify
% the Clustering Structure"

% Inputs:
% order: object ordering produced by OPTICS
% RD: object reachability distance
% CD: objct core distance
% eps: any eps less than generating eps
% geneps: generating eps value

% Outputs:
% Y: Cluster assignment

function [Y1,T,CDbw] = ExtractOpticsPartition(X1,X2,RD,CD,order,eps,geneps)
    if eps > geneps
        display('error');
    else
        noise = -1;
        clusterid = 0;
        for i = 1:length(order)
            if RD(i) > eps
                % if reachability distance is larger than eps, cluster ID is undefined 
                if CD(i) <= eps
                    clusterid = clusterid+1;
                    Y(i) = clusterid;
                else
                    Y(i) = noise;
                end
            else
                Y(i) = clusterid;
            end
        end     
    end
    
[~,ClusterSize] = Size_Analysis(Y);
percent = 2.0; % threshold
RC = ClusterSize(ClusterSize(:,3)>percent,1);
RC = sort(RC);
if ~isempty(Y==-1) && isempty(RC == -1)
    RC = [-1;RC];
end
Y1 = -1*ones(length(X1),1);

%Renames clusters with sequential numbers now that clusters with very few
%points in them have been reassigned to noise
if any(RC == -1) == 1
    Y1(Y==-1) = -1;
    for i = 2:length(RC)
        Y1(Y==RC(i)) = i-1;
    end
else
    for i = 1:length(RC)
        Y1(Y==RC(i)) = i;
    end
end

%Re-assign "low intensity clusters" (those consisting entirely of points
%with counts <6) to noise cluster
Y1 = LowI(X1(order,:),Y1);

[T,~] = Size_Analysis(Y1);

CDbw = CDbwIndex(X2(order,:),Y1);
end