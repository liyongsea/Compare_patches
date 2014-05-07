function []=showPatches(patches,nb)
% show the nb*nb first components of patches
% patches should be a (a*a)*d matrix
% where a*a is the size of a patch
    it=1;
    p_size=sqrt(size(patches,1));
    if nargin<2
        nb=p_size;
    end
    for i=1:nb*nb
        subplot(nb,nb,it),imshow(reshape(patches(:,i),p_size,p_size),[min(patches(:,i)),max(patches(:,i))]);
        it=it+1;
    end
end