%% read data
close all
clear all
rng('default');
%%
dataPath='../../data/%s';
%I=imread('/home/li/MVA/Graphcut_shadow/data/toulouse1_qb.gif');
% I=imread(sprintf(dataPath,'zebra.jpg'));
I= load_image('lena');
% I=rgb2gray(I);
%I=imresize(I,0.5);
% I=double(I)/255;
% I=0.5+randn(size(I));
n = 128;
c = [100 200];
I = rescale( crop(I,n, c) );
hf=figure(),imshow(I);
% set(gca,'position',[0 0 1 1],'units','normalized')
% saveas(hf,'wn','png');
%% extract patches
p_size=7; %odd
patches=extractPatches(I,p_size*p_size);
hf=figure,
[patches_principle,coeff_proj]=princomp(patches);
showPatches(patches_principle,4);
saveas(21,'patches','png');
%% distribution of the first components
hf=figure,
it=1;
for i=35:40
subplot(2,3,it) 
hist(coeff_proj(:,i),255)
it=it+1;
end
%set(gca,'position',[0 0 1 1],'units','normalized')
% saveas(hf,'hist','png');
%% approximate the hist distribution by a gaussian
pi=35;
Y=coeff_proj(:,pi);
Z=patches_principle(:,pi);
m=estimateGaussian(Y);
l=estimateLaplace(Y);
x=[-0.5:0.01:0.5]';
hf=figure,[h1,c1]=hist(Y,255)
hold on,plot(c1,h1./trapz(c1,h1),'-b','LineWidth',2)
saveas(hf,sprintf('h_%d',pi),'png');
% m.mu=0.8;
% drawGaussianMixture(m,x);
% m.mu=-0.2;
drawGaussianMixture(m,x);
drawLaplace(l,x);
saveas(hf,sprintf('h_%d_gau',pi),'png');
hf=figure,
imshow(reshape(Z,p_size,p_size),[min(Z),max(Z)]);
%set(gcf,'position',get(0,'screensize'))
set(gca,'position',[0 0 1 1],'units','normalized')
saveas(hf,sprintf('patch_%d',pi),'png');
%% estimate componant's distibution
resolution = 100;
[h c]=estimCom(coeff_proj,resolution);
[mu var]=estimGaussianCom(coeff_proj);
%% compare patches
P=patches_principle';
x1=zeros(p_size*p_size,1);
x1(1)=1;
x2=zeros(p_size*p_size,1);
x2(2)=1;
compareproba( x1,x2,coeff_proj,0.03,255 )
%%
% a contrario pior
priorModel.h=h;
priorModel.c=c;
% pior approximated by a gaussian
priorModel_gau.mu=mu;
priorModel_gau.var=var;
para.sigma=0.5;
px1=patches_principle*x1;
px2=patches_principle*x2;
similarity(px1,px2,P,priorModel,para)
similarity_gau_approx(px1,px2,P,priorModel_gau,para)
%%