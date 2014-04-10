function [ A ] = show_best_matches( I, p1, compare )
%SHOW_BEST_MATCHES shows best patches for a specific input method
%   arguments :
%       I - original input image
%       p1 - reference patch - should be square and side size should be odd
%       compare - function to compare the patches
%   return :
%       A - the matrix that contains all comparison results

figure,imshow(p1);
title('original patch');

p_size = size(p1,1);
ds=(p_size-1)/2;
A=zeros(size(I,1)-p_size+1,size(I,2)-p_size+1);
for i=1+ds:size(I,1)-ds
%    i
    for j=1+ds:size(I,2)-ds
        pij=I(i-ds:i+ds,j-ds:j+ds);
        A(i-ds,j-ds)=compare(p1(:),pij(:));
    end
end

figure,imshow(A,[min(A(:)), max(A(:))]);
axis on
colormap hot;
title('matches map');

[~,u]=sort(-A(:));
us=u(1:100);
[ys,xs]=ind2sub(size(A),us);
figure,imshow(I(1+ds:end-ds,1+ds:end-ds));
hold on
scatter(xs,ys,7,'r');
title('best matches');

figure,
title('best paches');
for i=1:16
    subplot(4,4,i);
    [y,x]=ind2sub(size(A),u(i));
    imshow(I(y:y+2*ds,x:x+2*ds));
end

end

