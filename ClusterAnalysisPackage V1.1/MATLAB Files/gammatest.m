function [normalizedgamma,gammaMC] = gammatest(Y1,X1,n_MC)

h = waitbar(0,'Please wait...');

n = length(Y1);
m = n*(n-1)/2;

% Creating assignment matrix
Y = zeros(n);
for i = 1:n-1
    for j = i+1:n
        if Y1(i) == Y1(j)
            Y(i,j) = 0;
        else
            Y(i,j) = 1;
        end
    end
end

% calculating proximity matrix
D = pdist(X1);
X = squareform(D);
X = triu(X);
%X = pdist2(X1,X1);
%X = triu(X);

% calculating gamma statistic
mx = 1/m*sum(sum(X));
my = 1/m*sum(sum(Y));
sx = sqrt(1/m*sum(sum(X.^2)) - mx^2);
sy = sqrt(1/m*sum(sum(Y.^2)) - my^2);

%gamma = sum(sum(triu(X.*Y,1)));
normalizedgamma = 1/m*(sum(sum(triu((X-mx).*(Y-my),1))))/(sx*sy);

waitbar(1/(n_MC+1))

for k = 1:n_MC

    waitbar((k+1)/(n_MC+1))

    % permuting cluster membership
    idx = randperm(length(Y1));

    % creating new assignment matrix
    for i = 1:n-1
        for j = i+1:n
            if Y1(idx(i)) == Y1(idx(j))
                Y(i,j) = 0;
            else
                Y(i,j) = 1;
            end
        end
    end

    gammaMC(k) = 1/m*(sum(sum(triu((X-mx).*(Y-my),1))))/(sx*sy);
end
close(h)

figure()
hist(gammaMC)
xlabel('Value of Gamma Statistic','fontsize',16) % x-axis label
ylabel('Frequency','fontsize',16) % y-axis label
title('Gamma Statistic Sampling Distribution (from MC)','fontweight','bold','fontsize',24)

end
