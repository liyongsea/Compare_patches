close all
clear all
%% patches in images
n = 128;
c = [100 200];
f0 = load_image('lena');
f0 = rescale( crop(f0,n, c) );
hf=figure(101);
imageplot(f0);

%% add noise
sigma = 0.1;
f = f0 + randn(n,n)*sigma;
hf=figure(202);
imageplot(clamp(f));

%% PCA
options.k = 3;          % half size for the windows
options.T = 0.06;       % width of the gaussian, relative to max(M(:))  (=1 here)
options.max_dist = 15;  % search width, the smaller the faster the algorithm will be
options.ndims = 49;     % number of dimension used for distance computation (PCA dim.reduc. to speed up)
options.do_patchwise = 1;


%%
M=f;
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
%% prior model 
H0 = perform_lowdim_embedding(f0,options);
w=2*options.k+1;
resh=@(h) reshape(h,[],w*w);
va=var(resh(H0));
mu=mean(resh(H0));
gamma=(sigma^2)./max(va,1e-9);
gamma(1)=0;
C=1./(2*sigma^2*(1+gamma).*(2+gamma));

%% assemble parameters
gauModel=[mu;gamma;C]';
[M1,Wx,Wy] = my_nlmeans_mex(Ma,H,Ha,Vx-1,Vy-1,T,max_dist, do_median, do_patchwise, mask_process, mask_copy, exclude_self,gauModel);


