close all
clear all
I=imread('../../data/lena.bmp');
I=double(rgb2gray(I))/255.;
%% noise image
%I=I+randn(size(I))/10;
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
para.sigma=0.01;
%% extract an original patch
p1=I(323:329,323:329);
figure,imshow(p1);
%% show best matches found in reference image
compare_sim=@(p,q)(similarity(p(:),q(:),patches_principle',priorModel,para));
% Aq=show_best_matches(I,p1,compare_sim);
% 
compare_euclid=@(p,q)(-sqrt(sum(sum((p-q).^2))));
% As=show_best_matches(I,p1,compare_euclid);
%% Performance on PFA
dico=load('../../data/dico_for_gamma_L1.62.mat');
dico=dico.dico;
dico=dico/255;
nb=196;
%% Performance 
sigma=0.3;
dico1=dico+randn(size(dico))*sigma;
dico1(dico1<0)=0;
dico1(dico1>1)=1;
dico2=dico+randn(size(dico))*sigma;
dico2(dico2<0)=0;
dico2(dico2>1)=1;
M=[];
figure,
count=1;
sigma_list=[0.05,0.07,0.1,0.15,0.2,0.3];
for ps=sigma_list
    para.sigma = ps;
    compare_sim=@(p,q)(similarity(p(:),q(:),patches_principle',priorModel,para));
   
    hold on,
    M(:,:,count)=dico_curves( dico1,dico2,compare_sim,[0 0 1] );
    count=count+1;
end
%%
hold on,
M_eu=dico_curves( dico1,dico2,compare_euclid,[1 0 0] );

%%
FA_sim=[];
TD_sim=[];
for i=1:length(sigma_list)
    [fa,td]=getPFA(M(:,:,i));
    FA_sim(:,i)=fa;
    TD_sim(:,i)=td;
end
[FA_eu,TD_eu]=getPFA(M_eu);
%%
figure,hold on,
plot(FA_eu,TD_eu,'-r');
plot(FA_sim(:,1),TD_sim(:,1),'-b','Linewidth',2);
plot(FA_sim(:,2),TD_sim(:,2),'-','Linewidth',2,'Color',[0 1 1]);
plot(FA_sim(:,3),TD_sim(:,3),'-','Linewidth',2,'Color',[1 1 0]);
plot(FA_sim(:,4),TD_sim(:,4),'-','Linewidth',2,'Color',[0 0 0]);
plot(FA_sim(:,5),TD_sim(:,5),'-','Linewidth',2,'Color',[0.1 0.1 0.1]);
plot(FA_sim(:,6),TD_sim(:,6),'-','Linewidth',2,'Color',[0.25 0.25 0]);
legend('euclidian','sigma = 0.05','sigma = 0.07','sigma = 0.1','sigma = 0.15','sigma = 0.2','sigma = 0.3');