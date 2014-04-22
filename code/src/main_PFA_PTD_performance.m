close all
clear all
rng('default')
%% read image
%I=imread('../../data/texture.jpg');
I=imread('../../data/zebra.jpg');
I=double(rgb2gray(I))/255.;
%I=0.5+randn(size(I));
figure,imshow(I)
%% learn by pca
p_size=8;
patches=extractPatches(I,p_size*p_size);
[patches_principle,coeff_proj]=princomp(patches);
showPatches(patches_principle);
%%
resolution=255;
[h c]=estimCom(coeff_proj,resolution);
[mu var]=estimGaussianCom(coeff_proj);
priorModel.h=h;
priorModel.c=c;
priorModel_gau.mu=mu;
priorModel_gau.var=var;

%% prepare noisy dictionary
dico=load('../../data/dico_for_gamma_L1.62.mat');
dico=dico.dico;
dico=dico/255;
nb=size(dico,1);

sigma=0.3;
dico1=dico+randn(size(dico))*sigma;
dico1(dico1<0)=0;
dico1(dico1>1)=1;
dico2=dico+randn(size(dico))*sigma;
dico2(dico2<0)=0;
dico2(dico2>1)=1;
%% compute the simalrity performance on different sigma
M=[];% M is the correlation matrix
figure,
count=1;
%sigma_list=[0.05,0.07,0.1,0.15,0.2,0.3];
sigma_list=[0.15];
for ps=sigma_list
    para.sigma = ps;
    %compare_sim=@(p,q)(similarity(p(:),q(:),patches_principle',priorModel,para));
    compare_sim=@(p,q)(similarity_gau_approx(p(:),q(:),patches_principle',priorModel_gau,para));
    hold on,
    M(:,:,count)=dico_curves( dico1,dico2,compare_sim,[0 0 1] );
    count=count+1
end
%%
hold on,
%compare_euclid=@(p,q)(-sqrt(sum(sum((p-q).^2))));
d=14;
compare_reduc=@(p,q)(-distance_reduc(p(:),q(:),d,patches_principle'));
M_eu=dico_curves( dico1,dico2,compare_reduc,[1 0 0] );

%%
FA_sim=[];
TD_sim=[];
for i=1:length(sigma_list)
    [fa,td]=getPFA(M(:,:,i));
    FA_sim(:,i)=fa;
    TD_sim(:,i)=td;
end
[FA_eu,TD_eu]=getPFA(M_eu);
%% plot
hf=figure,hold on,
plot(FA_eu,TD_eu,'-r','Linewidth',2);
plot(FA_sim(:,1),TD_sim(:,1),'-b','Linewidth',2);
% plot(FA_sim(:,2),TD_sim(:,2),'-','Linewidth',2,'Color',[0 1 1]);
% plot(FA_sim(:,3),TD_sim(:,3),'-','Linewidth',2,'Color',[1 1 0]);
% plot(FA_sim(:,4),TD_sim(:,4),'-','Linewidth',2,'Color',[0 0 0]);
% plot(FA_sim(:,5),TD_sim(:,5),'-','Linewidth',2,'Color',[0.1 0.1 0.1]);
% plot(FA_sim(:,6),TD_sim(:,6),'-','Linewidth',2,'Color',[0.25 0.25 0]);
% legend('euclidian','sigma = 0.05','sigma = 0.07','sigma = 0.1','sigma = 0.15','sigma = 0.2','sigma = 0.3');
saveas(hf,'perf','png');