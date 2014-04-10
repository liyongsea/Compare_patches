I=imread('../../data/lena.bmp');
I=double(rgb2gray(I))/255.;
%%
I=I+randn(size(I))/10;
I(I<0)=0; I(I>1)=1;
%%
figure,imshow(I);

%% extract patches
p_size=7;
patches=extractPatches(I,p_size*p_size);

[patches_principle,coeff_proj]=princomp(patches);
showPatches(patches_principle);
%%
resolution=255;
[h c]=estimCom(coeff_proj,resolution);
%%
priorModel.h=h;
priorModel.c=c;
para.sigma=0.01;
%%
p1=I(323:329,323:329);
figure,imshow(p1);
%%
compare=@(p,q)(similarity(p(:),q(:),patches_principle',priorModel,para));
Aq=show_best_matches(I,p1,compare);
%%
compare=@(p,q)(-sqrt(sum(sum((p-q).^2))));
As=show_best_matches(I,p1,compare);
%%
dico=load('../../data/dico_for_gamma_L1.62.mat');
dico=dico.dico;
dico=dico/255;
nb=196;
%%
dico1=dico+randn(size(dico))/3;
dico1(dico1<0)=0;
dico1(dico1>1)=1;
dico2=dico+randn(size(dico))/3;
dico2(dico2<0)=0;
dico2(dico2>1)=1;
%%
M=zeros(nb,nb);
for i=1:nb
    i
    for j=1:nb
        %M(i,j)=similarity(dico1(:,i),dico2(:,j),patches_principle',priorModel,para);
        M(i,j)=-sum((dico1(:,i)-dico2(:,j)).^2);
    end
end
%%
[~,idxs]=sort(M(:));
Id=eye(nb);
BI=Id(idxs);
TD=(1:numel(idxs))'-cumsum(BI);
FA=cumsum(BI);
figure,plot(FA/nb,TD/(nb*(nb-1)),'r');