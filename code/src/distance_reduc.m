function res=distance_reduc(x1,x2,d,P)
% compute the likelyhood ratio using the prior Model
%------------------
%projet x1 x2 (column vectors) on the principal componant basis
cx1=P*x1;
cx2=P*x2;

res=norm(cx1(1:d)-cx2(1:d));
end