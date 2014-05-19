function []=drawLaplace(l,x)
    y=1./(2*l.b)*exp(-abs(x-l.mu)/l.b);
    plot(x,y, '-g','LineWidth',1);
end