function [mu var]=estimGaussianCom(coeff_proj)
% estimate the probabilty distribution of each componant from their histogram
% h(i) : P(c=ci)
mu=[];
var=[];
    for i=1:size(coeff_proj,2)
        m=estimateGaussian(coeff_proj(:,i));
        mu=[mu,m.mu];
        var=[var,m.sigma]
    end
end