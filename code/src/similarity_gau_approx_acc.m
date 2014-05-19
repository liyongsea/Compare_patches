function res=similarity_gau_approx_acc(x1,x2,priorModel_gau_appro,para)
% compute the likelyhood ratio using the prior Model
% warning this is a accelarate version, some parameters must be precomputed
% and x1 x2 are already projected on the principale component basis
%------------------

% var_theta=reshape(priorModel_gau_appro.var',size(x1));
% mu_theta=reshape(priorModel_gau_appro.mu',size(x1));
% a=reshape(var_theta./((para.sigma)^2),size(x1));
var_theta=priorModel_gau_appro.vt;
mu_theta=priorModel_gau_appro.mt;
a=priorModel_gau_appro.a;

theta12=(a.*(x1+x2)+mu_theta)./(1+2*a);
theta1=(a.*x1+mu_theta)./(1+a);
theta2=(a.*x2+mu_theta)./(1+a);
num=(-0.5/para.sigma^2 * (x1-theta12).^2-0.5/para.sigma^2 * (x2-theta12).^2-0.5./var_theta.*(theta12-mu_theta).^2);
denom1=(-0.5/para.sigma^2 * (x1-theta1).^2-0.5./var_theta.*(theta1-mu_theta).^2);
denom2=(-0.5/para.sigma^2 * (x2-theta2).^2-0.5./var_theta.*(theta2-mu_theta).^2);

res=(num)-(denom1)-(denom2);
end