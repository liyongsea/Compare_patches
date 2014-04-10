function res=similarity_gau_approx(x1,x2,P,priorModel_gau_appro,para)
% compute the likelyhood ratio using the prior Model
%------------------
%projet x1 x2 (column vectors) on the principal componant basis
cx1=P*x1;
cx2=P*x2;
var_theta=priorModel_gau_appro.var';
mu_theta=priorModel_gau_appro.mu';
a=var_theta./((para.sigma)^2);
theta12=(a.*(cx1+cx2)+mu_theta)./(1+2*a);
theta1=(a.*cx1+mu_theta)./(1+a);
theta2=(a.*cx2+mu_theta)./(1+a);
num=sum(-0.5/para.sigma^2 * (cx1-theta12).^2-0.5/para.sigma^2 * (cx2-theta12).^2-0.5./var_theta.*(theta12-mu_theta).^2);
denom1=sum(-0.5/para.sigma^2 * (cx1-theta1).^2-0.5./var_theta.*(theta1-mu_theta).^2);
denom2=sum(-0.5/para.sigma^2 * (cx2-theta2).^2-0.5./var_theta.*(theta2-mu_theta).^2);
res=sum(num(:))-sum(denom1(:))-sum(denom2(:));
end