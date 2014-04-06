function [ patches ] = extractPatches( im,s )
%EXTRACTPATCHES Summary of this function goes here
%   Detailed explanation goes here
l=sqrt(s);
u=(l-1)/2;
x=size(im,1);
y=size(im,2);
z=1;
patches=zeros(s,(x-2*u)*(y-2*u));
for i=1+u:x-u
    for j=1+u:y-u
        p=im(i-u:i+u,j-u:j+u);
        patches(:,z)=p(:);
        z=z+1;
    end
end
patches=patches';
end

