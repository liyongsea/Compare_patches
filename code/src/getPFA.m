function [FA,TD]=getPFA(M)
nb=size(M,1);
[~,idxs]=sort(M(:));
Id=eye(nb);
BI=Id(idxs);
TD=(1:numel(idxs))'-cumsum(BI);
FA=cumsum(BI);
FA=FA/nb;
TD=TD/(nb*(nb-1))
end