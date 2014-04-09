I=imread('../../data/lena.bmp');
I=double(rgb2gray(I))/255.;
figure,imshow(I);
%% extract patches
p_size=7; %odd
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
figure,
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
figure,imshow(I);
hold on
scatter(xs,ys,7,'r');
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