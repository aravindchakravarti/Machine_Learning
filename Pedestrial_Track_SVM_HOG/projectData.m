function Z_data = projectData(X_data, U_data, K_vectors)
% Basically it is doing DOT PRODUCT. 
% This function Credit goes to Andrew Ng's Machine Learning team of
% Coursera

Z_data = zeros(size(X_data, 1), K_vectors);

first_k_columns_U = U_data(:,1:K_vectors);

Z_data = X_data * first_k_columns_U;

end