alist=0.0001:0.001:10;
f
a=0.5;
x_offl=-0.1:0.01:0.1
x1=0.5;
a_list=[0.1:0.5:10]
hf=figure
for a=a_list
d=[];
fa=[];
for diff=x_offl

x2=x1+diff;
mu=0;


theta=(a*(x1+x2)+mu)/(1+2*a);
theta1=(a*(x1)+mu)/(1+a);
theta2=(a*(x2)+mu)/(1+a);
f=(x1-theta)^2+(x2-theta)^2+1/a*(theta-mu)^2-(x1-theta1)^2-1/a*(theta1-mu)^2-(x2-theta2)^2--1/a*(theta2-mu)^2;

fa=[fa,f];
d=[d,0.5*(x1-x2)^2];
end

plot(x_offl,fa,'LineWidth',2)
hold on, plot(x_offl,d,'-r','LineWidth',2)
%legend('max a posteriori','euclidian')
%set(gca,'position',[0 0 1 1],'units','normalized')
end
grid, axis equal
saveas(hf,'curve','png');
