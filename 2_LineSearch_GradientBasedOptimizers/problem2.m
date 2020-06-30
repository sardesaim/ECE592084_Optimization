%problem 2
clc; close all;
%set range of x from -1 to 1 
x = linspace(-1,1,50);
y = x;
% create a meshgrid for plotting the function
[x1,x2] = meshgrid(x,y);
% define the given function for 3d plotting
f = (x2-x1).^4+12.*x1.*x2-x1+x2-3;
box on; mesh(x1,x2,f); %plot the function
figure;
% plot contours of the function
contour(x1,x2,f, 20)