classdef DataObjectClass
   properties
      value
      processed
      index
      RD
      CD
   end
   methods
       function obj = setCD(obj,nbhd,minpts)
           if size(nbhd.Neighbors) < minpts
               obj.CD = 0; % object is not a core point at any epsilon less than geneps
           else
               sortd = sort(nbhd.Distance); % sorting distances in the neighborhood in ascending order
               % choosing the smallest distance that would generate a neighborhood with minpts in it
               obj.CD = sortd(minpts);
           end
       end
   end
end