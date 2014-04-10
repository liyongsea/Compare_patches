close all
clear all
I=imread('../../data/lena.bmp');
I=double(rgb2gray(I))/255.;
%% noise image
I=I+randn(size(I))/10;
I(I<0)=0; I(I>1)=1;
%%
figure,imshow(I);

%% extract patches
p_size=8;
patches=extractPatches(I,p_size*p_size);

[patches_principle,coeff_proj]=princomp(patches);
showPatches(patches_principle);
%% parameters for similarity method
resolution=255;
[h c]=estimCom(coeff_proj,resolution);
priorModel.h=h;
priorModel.c=c;
para.sigma=0.1;
%% extract an original patch
p1=I(323:329,323:329);
figure,imshow(p1);
%% show best matches found in reference image
compare_sim=@(p,q)(similarity(p(:),q(:),patches_principle',priorModel,para));
% Aq=show_best_matches(I,p1,compare_sim);
% 
compare_euclid=@(p,q)(-sqrt(sum(sum((p-q).^2))));
% As=show_best_matches(I,p1,compare_euclid);
