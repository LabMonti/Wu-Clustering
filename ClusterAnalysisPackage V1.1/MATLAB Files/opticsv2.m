% Optics Algorithm based on the paper "OPTICS: Ordering Points To Identify
% the Clustering Structure"
% Distances calculated using manhattan metric
% current run time is O(n^2) but can be reduced to O(n) with implementation
% of indexing. I have heard that this can be difficult to execute in matlab

% Inputs:
% D: The input data set to be clustered
% eps: The generating (maximum) neighborhood radius 
% minpts: minimum amount of points required to form a cluster

% Outputs:
% RD: The reachability distance of each object
% CD: The core distance of each object
% order: The order in which each object was processed

% Undefined = 0


function [RD,CD,order]=opticsv2(D,eps,minpts)

m = size(D,1);
% tells us whether an object has been processed
OP = zeros(m,1);
% the order variable will tell us in which order the data is processed 
order = zeros(m,1);
% CD will store all the core distances
CD = zeros(m,1);
% RD will store all of the reachability distances
RD = zeros(m,1);
% stores object reachability distance
RD1 = zeros(m,1);
% Objects which are directly density-reach-able from a current core object 
% are inserted into the seed-list for further expansion
seeds = [];
k = 1; % order index

h = waitbar(0,'Please wait...');

% This is the main loop of the program. It checks to see if an object has
% been processed. If it has not been processed, we expand the cluster order
% around that object
for i = 1:m
    if isempty(find(order==i, 1)) == 1
       [RD,CD,order,k,RD1,OP]=ExpandClusterOrder(D,eps,minpts,OP,RD,CD,order,RD1,i,k,seeds,h);                
    end
end

close(h);

end

function [RD,CD,order,k,RD1,OP] = ExpandClusterOrder(D,eps,minpts,OP,RD,CD,order,RD1,i,k,seeds,h)

% this should give us the indices of the minpts neighbors of point D(i,:)
% implement region query
d = pdist2(D(i,:),D);
nbhd = find(d < eps);

order(k) = i; % object processed is true
OP(i) = 1;

m = size(D,1);

% calculating core distance of D(i,:)
% is this setting the core distance of everything in neighbors?
if size(nbhd) < minpts
    CD(k) = 0;
else
    % sorting distances in the neighborhood in ascending order
    sortd = sort(d(nbhd));
    % choosing the smallest distance that would generate a neighborhood
    % with minpts in it
    CD(k) = sortd(minpts);
end
    
% Establishing reachability distance
RD(k) = inf;
RD1(i) = inf;
    
%if core distance is not undefined
if CD(k) ~= 0
    % generate the seeds list
    [RD1,OP, seeds] = orderseedsupdate(D,nbhd,RD1,CD,OP,i,k,seeds);
    while ~isempty(seeds)
        % increment the order index
        k = k + 1;
        if mod(k,1000) == 0
            waitbar((k+0.0)/m,h);
        end
        % next object is object with smallest reachability distance
        % Sorting seeds with lowest reachability at the top
        % I might be able to innovate a bit here and come up with a
        % way of dealing with ties
        seeds = sortrows(seeds,2);        
        currentp = seeds(1); 
        d = pdist2(D(currentp,:),D);
        nbhd = find(d<eps); 
        OP(currentp) = 1; % current object is processed
        order(k) = currentp; % storing order in which we process the points
        RD(k) = seeds(1,2); % updating reachability distance of current processed point
        % Calculating core distance of current p
        if size(nbhd) < minpts
            CD(k) = 0;
        else
            sortd = sort(d(nbhd));
            CD(k) = sortd(minpts);
        end
            
        % if the object is a core object, we will find density reachable
        % objects and place them into seeds
        if CD(k) ~= 0 
            [RD1,OP,seeds] = orderseedsupdate(D,nbhd,RD1,CD,OP,currentp,k,seeds);  
        end
        seeds = seeds(2:size(seeds,1),:); 
    end
end
    k = k + 1;
    if mod(k,1000) == 0
        waitbar((k+0.0)/m,h);
    end
end

function [RD1,OP,seeds] = orderseedsupdate(D,nbhd,RD1,CD,OP,j,k,seeds)
    % the error is occuring here somewhere
    cdist = CD(k);
    q = D(j,:);
    
    for i = 1:length(nbhd)
        % We need to update core distance and reachability distance
        
        if OP(nbhd(i)) == 0 % if object has not been processed
            % manual calculation of Euclidean distance
            newrdist = max(cdist,sqrt((D(nbhd(i),1)-q(1))^2 + (D(nbhd(i),2)-q(2))^2 + (D(nbhd(i),3)-q(3))^2));
            if RD1(nbhd(i),1) == 0 
                % Inserting object into seeds list
                seeds = vertcat(seeds,[nbhd(i),newrdist]);
                RD1(nbhd(i),1) = newrdist; % storing reachability dist in point index
            elseif newrdist < RD1(nbhd(i))
                % we will check if nbhd(i) is already in seeds
                % we might not need to check if seeds is empty
                if size(seeds,1) == 0 
                    seeds = vertcat(seeds,[nbhd(i),newrdist]);                 
                else
                    % logical indexing for speed (suggested by matlab)
                    % we must replace the reachability distance here
                    seeds(seeds(:,1)==nbhd(i),2) = newrdist;
                    RD1(nbhd(i),1) = newrdist;
                end
            end
        end
    end 
end







