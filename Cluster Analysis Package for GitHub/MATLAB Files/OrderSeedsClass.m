classdef OrderSeedsClass
   properties
       List
       Dist
   end
   methods
       function [seed,S] = update(seed,S,nbhd,obj)
           c_dist = obj.CD;
           for j = 1:length(nbhd.Neighbors)
               object = nbhd.Neighbors(j);
               if object.processed == 0
                   new_r_dist = max(c_dist,sum(abs(object.value - obj.value)));
                   if isempty(object.RD) == 1 % if object reachability distance is not in seeds list
                       object.RD = new_r_dist; % setting object reachability distance
                       S.Set(object.index) = object; % updating database
                       
                       % inserting object into seedlist
                       seed.List = [seed.List;object];
                       seed.Dist = [seed.Dist;new_r_dist];
                   elseif new_r_dist < object.RD % if object is already in seeds list
                       object.RD = new_r_dist; % set object reachability distance to new_r_dist
                       S.Set(object.index) = object; % update database
                       
                       % updating object reachability distance in seed list
                       seedindex = [seed.List.index]; 
                       index = seedindex == object.index;
                       seed.List(index) = object;
                       seed.Dist(index) = new_r_dist;
                   end
               end
           end
       end
       
       function obj = next(seed) % chooses object in seed list with minimum reachability distance
           [~,idx] = sort(seed.Dist); 
           obj = seed.List(idx(1));
       end
   end
end