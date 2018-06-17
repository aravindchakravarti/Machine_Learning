function [U, S] = pca_manual(X)

[m, n] = size(X);

U = zeros(n);
S = zeros(n);

Sigma = (1/m)*(X'*X);

[U, S, V] = svd(Sigma);

end
