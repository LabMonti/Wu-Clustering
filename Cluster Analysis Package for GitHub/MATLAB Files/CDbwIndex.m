% Density Based Clustering Validation

function [CDbw] = CDbwIndex(X2,Y)
[idx,cnames] = grp2idx(Y);
n = length(idx);
k = length(cnames);
mbrs = (repmat(1:k,n,1) == repmat(idx,1,k));
Xc = cell(k,1);
rpoint = cell(k,1);
crp = cell(k,k); % closest representative points
crd = cell(k,k); % closest representataive point distance
rcr = cell(k,k); % respective closest representative points
rcd = cell(k,k);
dens = zeros(k,k);
dist = zeros(k,k);


% Finding representative points for our cluster
for j = 1:k
    Xc_j = X2(mbrs(:,j),:); % find members in cluster j
    stdev(j) = norm(std(Xc_j));
    Xc{j} = Xc_j;
    nc = size(Xc_j,1); % compute length of cluster
    r = ceil(.3*nc); % determine number of representative points as 30% of cluster size
    rp = []; % initialize array for storing representative points
    d = []; % initialize distance array
    dmin = pdist2(Xc_j,mean(Xc_j,1));
    for i = 1:r
        [~,I] = max(dmin); % find index of point furthest away from representative point
        rp = [rp;Xc_j(I,:)]; % store representative point into rp array
        d(:,i) = pdist2(Xc_j,Xc_j(I,:)); 
        dmin = min(dmin,d(:,i));
    end
    rpoint{j} = rp;
end

% used to view representative points if desired
%{
fig = cell2mat(rpoint);
figure()
scatter(fig(:,1),fig(:,2));
%scatter3(fig(:,1),fig(:,2),fig(:,3));
view(2)
%}

%beginning validation calculations
% first determine the closest representative points
for i=1:k 
    for j=1:k
        if i < j
            closest_rep = [];
            closest_rep_dist = [];
            for l=1:size(rpoint{j},1)
                rep_dist=pdist2(rpoint{j}(l,:),rpoint{i});
                [d,I] = min(rep_dist,[],2);  % find index and distance of closest representatives in cluster j to a representative point in cluster i
                closest_rep = [closest_rep;[I,l]];
                closest_rep_dist = [closest_rep_dist;d];
            end
            crp{i,j} = closest_rep;
            crd{i,j} = closest_rep_dist;
        elseif i > j
            closest_rep = [];
            closest_rep_dist = [];
            for l=1:size(rpoint{j},1)
                rep_dist=pdist2(rpoint{j}(l,:),rpoint{i});
                [d,I] = min(rep_dist,[],2);  % find index of closest representative in cluster j to a representative point in cluster i
                closest_rep = [closest_rep;[l,I]]; % ensuring order of indices match when taking intersection below
                closest_rep_dist = [closest_rep_dist;d];
            end
            crp{i,j} = closest_rep;
            crd{i,j} = closest_rep_dist;
        end
    end
end



% finds the respective closest representative points by taking the
% intersection between the crp{i,j} and crp{j,i}
for i = 1:k
    for j = 1:k
        [rcr{i,j},ia,~] = intersect(crp{i,j},crp{j,i},'rows');
        rcd{i,j} = crd{i,j}(ia);
    end
end



% Once we have found the respective closest representative points, we will
% now compute the density between clusters
for i = 1:k
    for j = i+1:k % saves on computation
        n_rcr = size(rcr{i,j},1); % determining the number of closest representative pairs
        vi = rpoint{i}(rcr{i,j}(:,1),:); % extracting representative points for cluster i
        vj = rpoint{j}(rcr{i,j}(:,2),:); % extracting representative points for cluster j
        u = .5*(vi+vj); % finds the midpoint between pairs of respective closest representative points
        c_ij = vertcat(Xc{i},Xc{j}); % storing points for clusters i and j
        for l = 1:size(u,1) % for each midpoint, ul
            d1 = pdist2(c_ij,u(l,:)); % calculate distance between a point in a cluster and midpoint ul
            f = d1 < (stdev(i)+stdev(j))/2; % finds points that are less than std away
            card = sum(f)/size(c_ij,1); % compute cardinality
            densl(l) = rcd{i,j}(l)/(stdev(i)+stdev(j))*card; % compute density
        end
        dist(i,j) = sum(rcd{i,j})/n_rcr; % average of all representative point distances
        dens(i,j) = sum(densl)/n_rcr;
        clear densl;
    end
end

% The intercluster separation computation
dens = dens + transpose(dens); % saves on computation, we can do this because being respective representative points is a symmetric relation
inter_dens = sum(max(dens,[],2))/k;
dist = dist + transpose(dist);
dist(dist==0) = inf; % so we do not choose minimum of 0
if size(dist) == 1
    sep = 0;
else
    sep = sum(min(dist,[],2))/k/(1+inter_dens);
end

% intracluster cohesion computation
stdev_a = mean(stdev);
for i = 1:8 % we will keep s, the shrinking factor, within [.1,.8] to avoid the trivial cases
    s = i*.1; % incrementing shrinking factor by 10% each iteration
    for j = 1:k
        v1 = rpoint{j}; % find the representative points for cluster j
        Xj = Xc{j}; % find all the points in cluster j
        center = mean(Xj); % find mean of cluster j and call this the center
        v = v1 + s* bsxfun(@minus,center(ones(size(v1,1),1),:),v1); % calculate the shrunken representative points
        d = pdist2(Xj, v); % calculate the pairwise distance between all points in cluster j and the shrunken representative points
        f = d < stdev_a; % find all points that are less than average std away from the shrunken representative points
        card = sum(f,1)/size(Xj,1); % compute cardinality of each shrunken representative point
        dens_cl(j) = sum(card)/size(v1,1); % sum up total and normalize
    end
    intra_dens(i) = sum(dens_cl)/k/stdev_a; % computer intra cluster density for different shrinking factors
end
compactness = sum(intra_dens)/8; % average the intra cluster densities

% we will now compute a metric for cluster cohesion. This is a measure for
% the density variation within a cluster
for i = 1:7
    intra_change(i) = intra_dens(i+1)-intra_dens(i); % computer difference in intracluster density 
end
intra_change = sum(intra_change)/7; % normalizing intra cluster density changes 

cohesion = compactness/(1+intra_change);
SC = sep*compactness;
CDbw = SC*cohesion;
end

