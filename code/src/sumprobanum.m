function [ num ] = sumprobanum( x1,x2,s,sigma2,res )
%SUMPROBANUM Summary of this function goes here
%   Detailed explanation goes here
num=0;
for i=1:numel(x1)
    [h,b]=hist(s(:,i),res);
    h=h/sum(h);
    num = num + log(exp((-1/sigma2)*(x1(i)-b).^2).*exp((-1/sigma2)*(x2(i)-b).^2)*h');
end
end

