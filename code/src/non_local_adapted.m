close all
clear all

%% patches in images

experience='../../results/byzh';%set path to save images
experience_noise='noise_05';
n = 128;
c = [120 200];
f0=double(rgb2gray(imread('../../data/byzh.jpg')))./255;
f0=imresize(f0,0.5);
% f0 = load_image('lena');
f0 = rescale( crop(f0,n, c) );
hf=figure(101);
imageplot(f0);
% set(gca,'position',[0 0 1 1],'units','normalized')
% saveas(hf,sprintf('%s/origine',experience),'png');

%% add noise
sigma = 0.1;
f = f0 + randn(n,n)*sigma;
hf=figure(202);
imageplot(clamp(f));
% set(gca,'position',[0 0 1 1],'units','normalized')
% saveas(hf,sprintf('%s/%s/noisy',experience,experience_noise),'png');

%% 
figure(),
w = 3; %half width
w1 = 2*w+1 ;

% location of pixels
[Y,X] = meshgrid(1:n,1:n);
% offsets
[dY,dX] = meshgrid(-w:w,-w:w);
% location of pixels to extract
dX = reshape(dX, [1 1 w1 w1]);
dY = reshape(dY, [1 1 w1 w1]);
X = repmat(X, [1 1 w1 w1]) + repmat(dX, [n n 1 1]);
Y = repmat(Y, [1 1 w1 w1]) + repmat(dY, [n n 1 1]);

X(X<1) = 2-X(X<1); Y(Y<1) = 2-Y(Y<1);
X(X>n) = 2*n-X(X>n); Y(Y>n) = 2*n-Y(Y>n);

patch = @(f)f(X + (Y-1)*n);%Patch extractor operator

P = patch(f);

clf;
for i=1:16
    x = floor( rand*(n-1)+1 );
    y = floor( rand*(n-1)+1 );
    imageplot( squeeze(P(x,y,:,:)), '', 4,4,i );
end
%% learn for L1
resh = @(P)reshape(P, [n*n w1*w1])';
[patches_principle,coeff_proj]=princomp(resh(patch(f0))');
figure,showPatches(patches_principle);
resolution=255;
priorModel_lap=estimateLaplace(coeff_proj);
para1.sigma=1*sigma;
para2.sigma=1*sigma;
%% learn for L2
[patches_principle,coeff_proj]=princomp(resh(patch(f0))');
figure,showPatches(patches_principle);
resolution=255;
[mu var]=estimGaussianCom(coeff_proj);
priorModel_gau.mu=mu;
priorModel_gau.var=var;
%%
iresh = @(P)reshape(P, [n,n w1*w1]);
[patches_principle,coeff_proj]=princomp(resh(patch(f))');
H=iresh(coeff_proj);
a=2;
H1=H(:,:,1:a);
H2=H(:,:,a+1:end);
%% L1 para
vresh =@(m) reshape(m,1,1,[]);
tronc1= @(m) m(:,:,1:a);
tronc2= @(m) m(:,:,a+1:end);
priorModel_lap.b=tronc1(vresh(priorModel_lap.b));
priorModel_lap.mu=tronc1(vresh(priorModel_lap.mu));
priorModel_lap.lambda=(para1.sigma^2)./priorModel_lap.b;
priorModel_lap.lambda(1)=0;
h=@(c,lambda)0.5*c.^2.*(abs(c)<=lambda)+(abs(c)>lambda).*(0.5*lambda.^2+lambda.*(abs(c)-lambda));
%% L2 para

priorModel_gau.vt=tronc2(vresh(priorModel_gau.var));
priorModel_gau.mt=tronc2(vresh(priorModel_gau.mu));
priorModel_gau.gamma=(para2.sigma^2)./priorModel_gau.vt;

priorModel_gau.C=1./(2*para2.sigma^2*(1+priorModel_gau.gamma).*(2+priorModel_gau.gamma));
%% mesure the similarity of one patch on all image
% L1 dist
repm=@(b)repmat(b,[n,n,1]);
mu=repm(priorModel_lap.mu);
lambda=repm(priorModel_lap.lambda);
compare_sim1=@(p,q)1/para1.sigma^2*sum(0.25*(p-q).^2+2*h((p+q)/2-mu,lambda/2)-h(p-mu,lambda)-h(q-mu,lambda),3);
distance1 = @(i)compare_sim1(H1,repm(H1(i(1),i(2),:)));
% L2 dist
gamma=repm(priorModel_gau.gamma);
C=repm(priorModel_gau.C);
m=repm(priorModel_gau.mt);
compare_sim2=@(p,q)sum(C.*((p-q).^2-2*gamma.*(p-m).*(q-m)),3);
distance2 = @(i)compare_sim2(H2,repm(H2(i(1),i(2),:)));

distance= @(i) distance1(i)+distance2(i);
normalize = @(K)K/sum(K(:));
kernel = @(i,tau)normalize( exp( -distance(i)/(2*tau^2) ) );
tau = 3;
i = [83 72]-58;
% i = [95 70];
D = distance(i);
K = kernel(i,tau);
figure(4),
imageplot(D, 'D', 1,2,1);
imageplot(K, 'K', 1,2,2);
%% NL filter

q = 15;
% L1
repm=@(b)repmat(b,[2*q+1,2*q+1,1]);
mu=repm(priorModel_lap.mu);
lambda=repm(priorModel_lap.lambda);
compare_sim1=@(p,q)1/para1.sigma^2*sum(0.25*(p-q).^2+2*h((p+q)/2-mu,lambda/2)-h(p-mu,lambda)-h(q-mu,lambda),3);

selection = @(i){clamp(i(1)-q:i(1)+q, 1,n), clamp(i(2)-q:i(2)+q,1,n)};
distance1 = @(i,sel)compare_sim1(H1(sel{1},sel{2},:),repm(H1(i(1),i(2),:)));
%     distance = @(i,sel)sum( abs(H(sel{1},sel{2},:) - repmat(H(i(1),i(2),:), ...
%         [length(sel{1}) length(sel{2}) 1])), 3 )/(w1*w1);
% L2
gamma=repm(priorModel_gau.gamma);
C=repm(priorModel_gau.C);
m=repm(priorModel_gau.mt);
compare_sim2=@(p,q)sum(C.*((p-q).^2-2*gamma.*(p-m).*(q-m)),3);
distance2 = @(i,sel)compare_sim2(H2(sel{1},sel{2},:),repm(H2(i(1),i(2),:)));

distance = @(i)distance1(i,selection(i))+distance2(i,selection(i));
kernel = @(i,tau)normalize( exp( -distance(i)/(2*tau^2) )/tau );
D = distance(i);
K = kernel(i,tau);
clf; 
imageplot(D, 'D', 1,2,1);
imageplot(K, 'K', 1,2,2);

%%
P3=reshape(patch(f),n,n,w1*w1);
NLval_P = @(K,sel)sum(sum(repmat(K,[1,1,size(P3,3)]).*P3(sel{1},sel{2},:,:)));
NLval_P = @(i,tau)NLval_P( kernel(i,tau), selection(i) );
tau_list=[0.5:0.5:4];
it=1;
mySNR=[];
gna_list=[];
g_list=[];
for tau=tau_list
h=zeros(size(P));
    for i=1:n
        if (mod(i,10)==1)
            fprintf('tau = %f, %f persent finished\n\r',tau,i/n*100)
        end
        for j=1:n
            h(i,j,:)=NLval_P([i,j],tau);
        end
    end
    gna_list(:,:,it)=h(:,:,w+1,w+1);
    g_list(:,:,it)=aggregation(h,7);
    mySNR(it)=snr(f0,g_list(:,:,it));
    it=it+1;
end
figure(),plot(tau_list,mySNR,'LineWidth',3);
title('SNR_tau');
%%
mySNR_w=[];
for i=1:length(tau_list)
    mySNR_w(i)=snr(f0,gna_list(:,:,i));  
end

hold on,plot(tau_list,mySNR_w,'-r','LineWidth',3);
title('SNR_tau');
%%
[snr_min,ind_opt]=max(mySNR);
%tau_list=0.004:0.002:0.015;
%for ind=1:length(tau_list)
ind=ind_opt;
tau_opt=tau_list(ind);
g_opt=g_list(:,:,ind);
hf=figure(155),imageplot(g_opt)
title(sprintf('SNR %f',snr(f0,g_opt)))

%%
% 
% options.k = 3;          % half size for the windows
% options.T = 0.06;       % width of the gaussian, relative to max(M(:))  (=1 here)
% options.max_dist = 15;  % search width, the smaller the faster the algorithm will be
% options.ndims = 30;     % number of dimension used for distance computation (PCA dim.reduc. to speed up)
% options.do_patchwise = 1;

