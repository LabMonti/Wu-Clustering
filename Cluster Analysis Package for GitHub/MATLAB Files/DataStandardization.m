function X1 = DataStandardization(X,w)

%We expect lowest log(G) to be negative since we are in units of G0
if min(X(:,2))>0
    X(:,2) = -X(:,2);
end
X1 = X;
X1(:,3) = 1./X1(:,3);
X1(:,1) = zscore(X1(:,1));
X1(:,2) = w*zscore(X1(:,2));
X1(:,3) = zscore(X1(:,3));

end