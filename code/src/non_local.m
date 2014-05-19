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

P = patch(f0);

clf;
for i=1:16
    x = floor( rand*(n-1)+1 );
    y = floor( rand*(n-1)+1 );
    imageplot( squeeze(P(x,y,:,:)), '', 4,4,i );
end

%% dimension reduction
d = 30;
resh = @(P)reshape(P, [n*n w1*w1])';
remove_mean = @(Q)Q - repmat(mean(Q), [w1*w1 1]);
P1 = remove_mean(resh(P));
C = P1*P1';
[V,D] = eig(C); D = diag(D);
[D,I] = sort(D, 'descend'); V = V(:,I);
figure(55),
plot(D); axis('tight');
%%
figure(55),
for i=1:16
    imageplot( reshape(V(:,i),[w1 w1]), '', 4,4,i );
end
%%
iresh = @(Q)reshape(Q', [n n d]);
descriptor = @(f)iresh( V(:,1:d)' * remove_mean(resh(P)) );
H = descriptor(f);
%%
distance = @(i)sum( (H - repmat(H(i(1),i(2),:), [n n 1])).^2, 3 )/(w1*w1);
normalize = @(K)K/sum(K(:));
kernel = @(i,tau)normalize( exp( -distance(i)/(2*tau^2) )/tau );
tau = .08;
%i = [83 72];
i=[73 76];
D = distance(i);
K = kernel(i,tau);
figure(55),
imageplot(D, 'D', 1,2,1);
imageplot(K, 'K', 1,2,2);
%% NL filter
q = 15;
selection = @(i){clamp(i(1)-q:i(1)+q, 1,n), clamp(i(2)-q:i(2)+q,1,n)};
distance = @(i,sel)sum( (H(sel{1},sel{2},:) - repmat(H(i(1),i(2),:), ...
        [length(sel{1}) length(sel{2}) 1])).^2, 3 )/(w1*w1);
%     distance = @(i,sel)sum( abs(H(sel{1},sel{2},:) - repmat(H(i(1),i(2),:), ...
%         [length(sel{1}) length(sel{2}) 1])), 3 )/(w1*w1);
distance = @(i)distance(i,selection(i));
kernel = @(i,tau)normalize( exp( -distance(i)/(2*tau^2) )/tau );
D = distance(i);
K = kernel(i,tau);
clf; 
imageplot(D, 'D', 1,2,1);
imageplot(K, 'K', 1,2,2);
%%
NLval = @(K,sel)sum(sum(K.*f(sel{1},sel{2})));
NLval = @(i,tau)NLval( kernel(i,tau), selection(i) );
[Y,X] = meshgrid(1:n,1:n);
NLmeans = @(tau)arrayfun(@(i1,i2)NLval([i1 i2],tau), X,Y);
h=NLmeans(0.011);
figure,imageplot(h)
%%

P3=reshape(patch(f),n,n,w1*w1);
NLval_P = @(K,sel)sum(sum(repmat(K,[1,1,size(P3,3)]).*P3(sel{1},sel{2},:,:)));
NLval_P = @(i,tau)NLval_P( kernel(i,tau), selection(i) );
tau_list=0.01;
it=1;
mySNR=[];
for tau=0.015
h=zeros(size(P));
    for i=1:n
        if (mod(i,10)==1)
            fprintf('tau = %f, %f persent finished\n\r',tau,i/128*100)
        end
        for j=1:n
            h(i,j,:)=NLval_P([i,j],tau);
        end
    end
    h_list(:,:,it)=aggregation(h,7);
    mySNR(it)=snr(f0,h_list(:,:,it));
    it=it+1;
end
figure(),plot(tau_list,mySNR,'LineWidth',3);
title('SNR_tau');
%%
h=[];
it=1;
tau_list=0.01:0.02:0.15;
mySNR=[];
for tau=tau_list
%tau = .06;
h(:,:,it)=NLmeans(tau);
mySNR(it)=snr(f0,h(:,:,it));
it=it+1;
end
figure(),plot(tau_list,mySNR,'LineWidth',3);
title('SNR_tau');

%%
[snr_min,ind_opt]=max(mySNR);
%tau_list=0.004:0.002:0.015;
%for ind=1:length(tau_list)
ind=2;
tau_opt=tau_list(ind_opt);
h_opt=h(:,:,ind_opt);
hf=figure(155),imshow(h_opt)
title(sprintf('SNR %f',snr(f0,h_opt)))
% set(gca,'position',[0 0 1 1],'units','normalized')
% saveas(hf,sprintf('%s/%s/eucli_reduc5_tau%g_SNR%g',experience,experience_noise,tau_opt,snr(f0,h_opt)),'png');
%end
%%
ind=6;
figure(),imageplot(h(:,:,ind))
title(sprintf('SNR %f',snr(f0,h(:,:,ind))))





%% my similarity for NL denoising
% learn by PCA
[patches_principle,coeff_proj]=princomp(resh(patch(f))');
figure,showPatches(patches_principle);
resolution=255;
[mu var]=estimGaussianCom(coeff_proj);
priorModel_gau.mu=mu;
priorModel_gau.var=var;
para.sigma=0.03;

%% assign value
priorModel_gau.vt=reshape(priorModel_gau.var,[],1);
priorModel_gau.mt=reshape(priorModel_gau.mu,[],1);
priorModel_gau.a=reshape(var./((para.sigma)^2),[],1);
compare_sim=@(p,q)(similarity_gau_approx_acc(p,q,priorModel_gau,para));

jresh=@(p)reshape(p,[n,n,w1*w1]);
J=jresh(coeff_proj);
J=permute(J,[3 1 2]);%J is permutated to fit the bsxfun (change it if you get a better idea)

%% mesure the similarity of one patch on all image

perm=@(k)permute(k,[2 3 1]);
distance = @(i)-perm(sum(bsxfun(compare_sim, J ,J(:,i(1),i(2)))));
normalize = @(K)K/sum(K(:));
kernel = @(i,tau)normalize( exp( -distance(i)/(2*tau^2) ) );
tau = 8;
%i = [83 72];
i = [95 70];
D = distance(i);
K = kernel(i,tau);
figure(4),
imageplot(D, 'D', 1,2,1);
imageplot(K, 'K', 1,2,2);

%% similarity kernel on a sub window

q = 14;
selection = @(i){clamp(i(1)-q:i(1)+q, 1,n), clamp(i(2)-q:i(2)+q,1,n)};
distance = @(i,sel)-perm(sum(bsxfun(compare_sim, J(:,sel{1},sel{2}),J(:,i(1),i(2)))));
distance = @(i)distance(i,selection(i));
kernel = @(i,tau)normalize( exp( -distance(i)/(2*tau^2) ) );
D = distance(i);
K = kernel(i,tau);
figure(4)
imageplot(D, 'D', 1,2,1);
imageplot(K, 'K', 1,2,2);

%% non local means denoising using our similarity criterion
NLval = @(K,sel)sum(sum(K.*f(sel{1},sel{2})));
NLval = @(i,tau)NLval( kernel(i,tau), selection(i) );
[Y,X] = meshgrid(1:n,1:n);
NLmeans = @(tau)arrayfun(@(i1,i2)NLval([i1 i2],tau), X,Y);
g=zeros(n,n);
g_list=[];
tau_list=5:0.2:7;
it=1;
mySNR=[];
for tau = tau_list
    for i=1:n
        if (mod(i,10)==1)
            fprintf('tau = %f, %f persent finished\n\r',tau,i/128*100)
        end
        for j=1:n
            g(i,j)=NLval([i,j],tau);

        end
    end
    mySNR(it)=snr(f0,g);
    g_list(:,:,it)=g;
    it=it+1;
end
figure(),plot(tau_list,mySNR,'LineWidth',3);
title('SNR_tau');
%% show the result with a best snr
[snr_min,ind_opt]=max(mySNR);
for ind=1:length(tau_list)
tau_opt=tau_list(ind);
g_opt=g_list(:,:,ind);
hf=figure(144),imshow(g_opt)
title(sprintf('SNR %f',snr(f0,g_opt)))
set(gca,'position',[0 0 1 1],'units','normalized')
saveas(hf,sprintf('%s/%s/our_tau%g_SNR%g',experience,experience_noise,tau_opt,snr(f0,g_opt)),'png');
end
%%
% 
% %%
% 
% 
% 
% 
% 
% 
% 
% 
% jresh=@(p)reshape(p,[n,n,w1*w1]);
% J=jresh(P);
% distance = @(i)-reshape( sum(bsxfun(compare_sim, resh(P), reshape(J(i(1),i(2),:),[w1*w1,1])  )),n,n);
% normalize = @(K)K/sum(K(:));
% kernel = @(i,tau)normalize( exp( -distance(i)/(2*tau^2) ) );
% tau = 7;
% i = [83 72];
% D = distance(i);
% K = kernel(i,tau);
% figure(4),
% imageplot(D, 'D', 1,2,1);
% imageplot(K, 'K', 1,2,2);
% %%
% 
% q = 14;
% selection = @(i){clamp(i(1)-q:i(1)+q, 1,n), clamp(i(2)-q:i(2)+q,1,n)};
% distance = @(i,sel)-reshape( sum(...
%     bsxfun(compare_sim, ...
%     reshape(P(sel{1},sel{2},:,:),[length(sel{1})*length(sel{2}),w1*w1])',...
%     reshape(J(i(1),i(2),:),[w1*w1,1])) )...
%     ,[length(sel{1}) length(sel{2})]);
% distance = @(i)distance(i,selection(i));
% kernel = @(i,tau)normalize( exp( -distance(i)/(2*tau^2) ) );
% D = distance(i);
% K = kernel(i,tau);
% figure(4)
% imageplot(D, 'D', 1,2,1);
% imageplot(K, 'K', 1,2,2);
% %%
% NLval = @(K,sel)sum(sum(K.*f(sel{1},sel{2})));
% NLval = @(i,tau)NLval( kernel(i,tau), selection(i) );
% [Y,X] = meshgrid(1:n,1:n);
% NLmeans = @(tau)arrayfun(@(i1,i2)NLval([i1 i2],tau), X,Y);
% g=zeros(n,n)
% tau = 7;
% for i=1:n
%     i
%     for j=1:n
%         g(i,j)=NLval([i,j],tau);
%         
%     end
% end
% figure(144)
% imageplot(g);