function res = evaluate_gaussian(data,distribution);
%% res = evaluate_gaussian(data,distribution);
%%
%% data: (dimension_data) X (number_data) array
%% distribution: structure with fields mn (mean) and cov (covariance)
%% 
%% returns res =  1/((2pi)^(n/2) det(cov)^(1/2)) exp( - 1/2 (x- mn)^T cov^{-1} (x-mn))
%%  
[ndims,ndata]  = size(data);
denominator    =  (2*pi)^(ndims/2)*sqrt(det(distribution.cov));
diff_from_mean = (data-repmat(distribution.mn,[1,ndata]));
precision_matrix =  inv(distribution.cov);

res  = 1/(denominator)*exp(-sum(diff_from_mean.*(precision_matrix*(diff_from_mean)),1)/2);