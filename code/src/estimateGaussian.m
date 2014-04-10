function m=estimateGaussian(x)
    m.mu=mean(x);
    m.sigma=var(x);
    m.pi=1;
end