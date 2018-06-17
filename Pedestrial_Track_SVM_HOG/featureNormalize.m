function [X_norm, mu_mean, sigma] = featureNormalize(X_data)

mu_mean = mean(X_data);
X_norm = bsxfun(@minus, X_data, mu_mean);

sigma = std(X_norm);
X_norm = bsxfun(@rdivide, X_norm, sigma);


end