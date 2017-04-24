classdef DataSetClass
   properties
       Set 
       Distance
       Neighbors
       X
   end
   methods
       function nbhd = getneighbors(S,obj,eps) % calculates epsilon nbhd of an object
           d = pdist2(obj.value,S.X,'cityblock');
           idx = find(d < eps & d > 0); % find all distinct points less than epsilon distance away 
           nbhd.Distance = d(idx);
           nbhd.Neighbors = S.Set(idx);
       end
   end
end