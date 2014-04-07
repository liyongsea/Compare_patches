%% read data
close all
clear all
rng('default');
%%
dataPath='../../data/%s';
%I=imread('/home/li/MVA/Graphcut_shadow/data/toulouse1_qb.gif');
I=imread(sprintf(dataPath,'zebra.jpg'));
I=rgb2gray(I);
I=imresize(I,0.5);
figure(),imshow(I);
I=double(I)/255;
%% extract patches
p_size=7; %odd
patches=extractPatches(I,p_size*p_size);

[patches_principle,coeff_proj]=princomp(patches);
showPatches(patches_principle);
%%
figure,hist(coeff_proj(:,1),255)
%% estimate componant's distibution
resolution = 255;
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