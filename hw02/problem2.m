%problem 2
clc; close all;
x = linspace(-1,1,50);
y = x;
[x1,x2] = meshgrid(x,y);
f = (x2-x1).^4+12.*x1.*x2-x1+x2-3;
box on; mesh(x1,x2,f);
figure;
contour(x1,x2,f, 20)