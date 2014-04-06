function [ r ] = compareproba( x1,x2,s,sigma2,res )
%COMPAREPROBA Summary of this function goes here
%   Detailed explanation goes here
r=sumprobanum(x1,x2,s,sigma2,res)-sumproba(x1,s,sigma2,res)-sumproba(x2,s,sigma2,res);

end

