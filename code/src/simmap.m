I=imread('../../data/lena.bmp');
I=double(rgb2gray(I))/255.;
%%
I=I+randn(size(I))/10;
I(I<0)=0; I(I>1)=1;
%%
figure,imshow(I);

%% extract patches
p_size=8;
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
ds=(p_size-1)/2;
As=zeros(size(I,1)-p_size+1,size(I,2)-p_size+1);
for i=1+ds:size(I,1)-ds
    i
    for j=1+ds:size(I,2)-ds
        pij=I(i-ds:i+ds,j-ds:j+ds);
        As(i-ds,j-ds)=similarity(p1(:),pij(:),patches_principle',priorModel,para);
    end
end
%%
figure,imshow(As,[min(As(:)), max(As(:))]);
axis on
%%
[~,u]=sort(-As(:));
s=16;
us=u(1:100);
[ys,xs]=ind2sub(size(As),us);
figure,imshow(I(1+ds:end-ds,1+ds:end-ds));
hold on
scatter(xs,ys,7,'r');
%%
figure,
for i=1:16
    subplot(4,4,i);
    [y,x]=ind2sub(size(As),u(i));
    imshow(I(y:y+2*ds,x:x+2*ds));
end
%%
figure,
ds=(p_size-1)/2;
Aq=zeros(size(I,1)-p_size+1,size(I,2)-p_size+1);
for i=1+ds:size(I,1)-ds
    i
    for j=1+ds:size(I,2)-ds
        pij=I(i-ds:i+ds,j-ds:j+ds);
        Aq(i-ds,j-ds)=-sqrt(sum(sum((pij-p1).^2)));
    end
end
figure,imshow(Aq,[min(Aq(:)), max(Aq(:))]);
%%
[~,u]=sort(-Aq(:));
s=16;
us=u(1:100);
[ys,xs]=ind2sub(size(Aq),us);
figure,imshow(I(1+ds:end-ds,1+ds:end-ds));
hold on
scatter(xs,ys,7,'r');
%%
figure,
for i=1:16
    subplot(4,4,i);
    [y,x]=ind2sub(size(Aq),u(i));
    imshow(I(y:y+2*ds,x:x+2*ds));
end
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