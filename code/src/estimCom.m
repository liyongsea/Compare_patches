function [h,c]=estimCom(coeff_proj,resolution)
% estimate the probabilty distribution of each componant from
% their histogram
% h(i) : P(c=ci)
h=zeros(resolution,size(coeff_proj,2));
c=zeros(resolution,size(coeff_proj,2));
    for i=1:size(coeff_proj,2)
        [hi,ci]=hist(coeff_proj(:,i),resolution);
        hi=hi./trapz(ci,hi);
        h(:,i)=hi;
        c(:,i)=ci;
    end
end