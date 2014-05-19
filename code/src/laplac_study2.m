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
sigma = 0.05;
% hf=figure(202);
% imageplot(clamp(f));
% set(gca,'position',[0 0 1 1],'units','normalized')
% saveas(hf,sprintf('%s/%s/noisy',experience,experience_noise),'png');

%% 
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


%% estimation on origin image
resh = @(P)reshape(P, [n*n w1*w1])';
[patches_principle,coeff_proj]=princomp(resh(patch(f0))');
% figure,showPatches(patches_principle);
resolution=255;
priorModel_lap0=estimateLaplace(coeff_proj);

para.sigma=0.4*sigma;
%% estimation on origin image

for simu=1:100
f = f0 + randn(n,n)*sigma;
[patches_principle,coeff_proj]=princomp(resh(patch(f))');
% figure,showPatches(patches_principle);
resolution=255;
% [mu var]=estimGaussianCom(coeff_proj);
% priorModel_lap.mu=mu;
% priorModel_lap.b=sqrt(max((var-sigma^2)/2,0));
priorModel_lap=estimateLaplace(coeff_proj);

%%
err(:,:,simu)=priorModel_lap0.b-priorModel_lap.b;
simu
end