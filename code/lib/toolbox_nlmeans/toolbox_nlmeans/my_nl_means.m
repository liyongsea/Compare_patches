function [M1,Wx,Wy,mask_copy] = my_nl_means(M, options)


options.null = 0;
Ma = getoptions(options, 'Ma', M);
if isempty(Ma)
    Ma = M;
end
T = getoptions(options, 'T', .05);
do_median = getoptions(options, 'do_median', 0);
do_patchwise = getoptions(options, 'do_patchwise', 0);
max_dist = getoptions(options, 'max_dist', 10);
mask_process = getoptions(options, 'mask_process', []);
mask_copy = getoptions(options, 'mask_copy', []);
exclude_self = getoptions(options, 'exclude_self', 0);

[m,n,s] = size(M);
[ma,na,sa] = size(Ma);

if isfield(options, 'Vx') && isfield(options, 'Vy')
    Vx = options.Vx;
    Vy = options.Vy;
else
    if na==n && ma==m
        [Vy,Vx] = meshgrid(1:n,1:m);
    else
        Vx = floor( rand(m,n)*(ma-1) ) + 1;
        Vy = floor( rand(m,n)*(na-1) ) + 1;
    end
end

% lift to high dimensional patch space
if not(isfield(options, 'Ha')) || not(isfield(options, 'P')) || not(isfield(options, 'Psi')) 
    [Ha,options.P, options.Psi] = perform_lowdim_embedding(Ma,options);
else
    Ha = options.Ha;
end
if not(isfield(options, 'H'))
    H = perform_lowdim_embedding(M,options);
else
    H = options.H;
end
[M1,Wx,Wy] = perform_nlmeans_mex(Ma,H,Ha,Vx-1,Vy-1,T,max_dist, do_median, do_patchwise, mask_process, mask_copy, exclude_self);
% convert back to matlab notation >0
Wx = Wx + 1;
Wy = Wy + 1;

if not(isempty(mask_process))
    I0 = find(mask_process<0.5); I = [];
    for j=1:s
        I = [I; I0+(j-1)*n*m];
    end
    M1(I) = M(I);
end
if not(isempty(mask_copy))
    mask_copy = 1-double(M1(:,:,1)~=-1);
    M1(M1==-1) = M(M1==-1);
end