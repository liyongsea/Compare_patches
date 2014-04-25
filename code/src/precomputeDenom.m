function D=precomputeDenom(J,priorModel_gau,para)
var_theta=priorModel_gau.vt;
mu_theta=priorModel_gau.mt;
a=priorModel_gau.a;
D=zeros(size(J));
for i=1:size(J,2)
    for j=1:size(J,3)
        x1=J(:,i,j);
        theta1=(a.*x1+mu_theta)./(1+a);
        D(:,i,j)=(-0.5/para.sigma^2 * (x1-theta1).^2-0.5./var_theta.*(theta1-mu_theta).^2);
    end
end
end