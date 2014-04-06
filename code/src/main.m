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