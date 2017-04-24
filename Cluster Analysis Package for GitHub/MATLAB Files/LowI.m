function Y=LowI(orderedX,Y)

[idx,cnames] = grp2idx(Y);
n = length(idx);
k = length(cnames);
mbrs = (repmat(1:k,n,1) == repmat(idx,1,k));

m = 1;
for i = 2:k
    if orderedX(mbrs(:,i),3) < 6
        Y(mbrs(:,i)) = -1;
    else
        Y(mbrs(:,i)) = m;
        m = m+1;
    end
end
end