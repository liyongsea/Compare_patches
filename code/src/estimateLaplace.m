function l=estimateLaplace(x)
    N=length(x);
    l.mu=median(x);
    l.b=1./N*sum(abs(x-repmat(l.mu,size(x,1),1)));
end