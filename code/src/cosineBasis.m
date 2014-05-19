function [B]=cosineBasis(n)
% get the n*n cosine basis
e=zeros(n,n);
B=zeros(n*n,n*n);
D=dctmtx(n);
for i=1:n*n
    k=e;
    k(i)=1;
    p=D*k*D';
    B(:,i)=p(:);
end
