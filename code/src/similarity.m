function res=similarity(x1,x2,P,priorModel,para)
% compute the likelyhood ratio using the prior Model
%------------------
%projet x1 x2 (column vectors) on the principal componant basis
cx1=P*x1;
cx2=P*x2;
diff_cx1=bsxfun(@minus,priorModel.c,cx1');
diff_cx2=bsxfun(@minus,priorModel.c,cx2');
exp_diff_cx1=exp(-1/(2*para.sigma*2)*(diff_cx1.^2));
exp_diff_cx2=exp(-1/(2*para.sigma*2)*(diff_cx2.^2));

denom=log(sum(exp_diff_cx1.*priorModel.h))+log(sum(exp_diff_cx2.*priorModel.h));
num=log(sum(exp_diff_cx1.*exp_diff_cx2.*priorModel.h));
res=sum(num(:))-sum(denom(:));
end