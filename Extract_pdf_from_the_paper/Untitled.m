clc;clear all;close all;

load carsmall

% plot(MPG);
x=0:100;
% sigma=1; mu=50;
% data= normpdf(x,mu,sigma);
% 
PD = fitdist(MPG, 'normal')

data= normpdf(x,23.7181,8.03573);
% 
% 
% P = normcdf(x,mu,sigma);
% 
% norminv(rand,mu,sigma)
% 
% 
% % [mu,sigma] = normfit(data)
% % 
% % Y = normpdf(x,mu,sigma);
% 
% figure(1);hold on
plot(x, data);
% plot(x,P);
% % plot(x,Y,'-r');
