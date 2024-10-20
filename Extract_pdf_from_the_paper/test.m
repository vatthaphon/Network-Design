function test
clc;clear all;close all

xdata = ...
 [0.9 1.5 13.8 19.8 24.1 28.2 35.2 60.3 74.6 81.3];
ydata = ...
 [455.2 428.6 124.1 67.3 43.2 28.1 13.1 -0.4 -1.3 -1.5];
x0 = [100; -1] % Starting guess
[x,resnorm] = lsqcurvefit(@myfun,x0,xdata,ydata)

F = x(1)*exp(x(2)*xdata);
figure(1);hold on
plot(xdata,ydata,'*r');
plot(xdata,F);
end

function F = myfun(x,xdata)
F = x(1)*exp(x(2)*xdata);
end