function [ M ] = dico_curves( dico1,dico2,compare,color )
%DICO_CURVES Summary of this function goes here
%   Detailed explanation goes here
nb=size(dico1,2);
M=zeros(nb,nb);
for i=1:nb
   i;
    for j=1:nb
        M(i,j)=compare(dico1(:,i),dico2(:,j));
    end
end
[~,idxs]=sort(M(:));
Id=eye(nb);
BI=Id(idxs);
TD=(1:numel(idxs))'-cumsum(BI);
FA=cumsum(BI);
plot(FA/nb,TD/(nb*(nb-1)),'Color',color);
end

