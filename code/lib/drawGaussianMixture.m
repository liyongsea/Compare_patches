function []=drawGaussianMixture(mixture,xSmp)
nb_ctrs=length(mixture.mu);
 [yMix yGau]=evaluate_mixture(xSmp',mixture);
hold on,
 for i=1:nb_ctrs
   plot(xSmp,yGau(:,i), '-b','LineWidth',1);
 end
  hold on,plot(xSmp,sum(yMix,2),'-r','LineWidth',2)
 
end