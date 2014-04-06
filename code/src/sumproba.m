function [ deno ] = sumproba( x,s,sigma2,res )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
deno=0;
for i=1:numel(x)
    [h,b]=hist(s(:,i),res);
    h=h/sum(h);
    deno = deno + log(exp((-1/sigma2)*(x(i)-b).^2)*h');
end
end

