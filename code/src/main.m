%% read data
close all
clear all
rng('default');
%%
dataPath='../../data/%s';
%I=imread('/home/li/MVA/Graphcut_shadow/data/toulouse1_qb.gif');
I=imread(sprintf(dataPath,'lena.bmp'));
I=rgb2gray(I);
%I=0.5+randn(size(I));
I=imresize(I,0.5);
figure(),imshow(I);
I=double(I)/255;
%% extract patches
p_size=7; %odd
patches=extractPatches(I,p_size*p_size);

[patches_principle,coeff_proj]=princomp(patches);
showPatches(patches_principle);
%%
figure,
for i=1:6
subplot(2,3,i) 
hist(coeff_proj(:,i)./sum(coeff_proj(:,i)),255)
end
%%
Y=coeff_proj(:,5);
m=estimateGaussian(Y);
x=[-2:0.01:2]';
figure,[h1,c1]=hist(Y,255)
hold on,plot(c1,h1./trapz(c1,h1),'-b')
drawGaussianMixture(m,x)
%% estimate componant's distibution
resolution = 100;
[h c]=estimCom(coeff_proj,resolution);

%% compare patches
P=patches_principle';
x1=zeros(p_size*p_size,1);
x1(1)=1;
x2=zeros(p_size*p_size,1);
x2(2)=1;
compareproba( x1,x2,coeff_proj,0.03,255 )
%%
priorModel.h=h;
priorModel.c=c;
para.sigma=0.5;
px1=patches_principle*x1;
px2=patches_principle*x2;
similarity(px1,px2,P,priorModel,para)