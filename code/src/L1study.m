clear all
close all
%%
lambda_list=0:1:15;
it=1;
for lambda=lambda_list
    
% lambda=0;
thresh=@(c,lambda)sign(c).*max(abs(c)-lambda,0);
h=@(c,lambda)0.5*c.^2.*(abs(c)<=lambda)+(abs(c)>lambda).*(0.5*lambda^2+lambda.*(abs(c)-lambda));
sq=@(c)0.5*c.^2;
c=[-5:0.01:5];
% figure(),plot(c,h(c,lambda))
hold on, plot(c,sq(c),'-r')

%%

x1=5;
x2=[-5:0.01:15];
dist=@(x1,x2,lambda)2*h((x1+x2)/2,lambda/2)-h(x1,lambda)-h(x2,lambda)+0.25*(x1-x2).^2;
figure(3),hold on, plot(x2,dist(x1,x2,lambda),'LineWidth',2,'Color',[it/length(lambda_list) 0 1-it/length(lambda_list)])
it=it+1;
end
hold on, plot(x2,0.25*(x1-x2).^2,'LineWidth',2,'Color',[1 0 0 ])
%%
ss=[6000,3000];
A=10*randn(ss);
B=10*randn(ss);
s=(A-B).^2;
ds=dist(A,B,lambda);