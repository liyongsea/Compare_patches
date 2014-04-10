function [yMix yGau]=evaluate_mixture(X,mixture)
 nb_ctrs=size(mixture.mu,1);
 yGau=zeros(size(X,2),nb_ctrs);
 for i=1:nb_ctrs
   distr.mn=mixture.mu(i,:)';
   distr.cov=mixture.sigma(:,:,i);
   yGau(:,i)=mixture.pi(i)*evaluate_gaussian(X,distr)';
 end
 yMix=max(yGau,[],2);
end