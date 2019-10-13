clc; close all;

%Gradient descent 
syms x1 x2;
f = (x2-x1).^4+12.*x1.*x2-x1+x2-3;
%Initializations for the sequence of points
X_init_array1{1} = [0.55, 0.7]; %starting point 1
X_init_array1{2} = [0,0];
X_init_array2{1} = [-0.9 -0.5]; %starting point 2
X_init_array2{2} = [0,0];
%Gradient descent from the two starting points 
X_init_array1 = gradDesc(f, X_init_array1);
X_init_array2 = gradDesc(f, X_init_array2);

%plotting the contours/ level sets for the function to plot sequence later
x = linspace(-1,1,50);
y = x;
[x1,x2] = meshgrid(x,y);
f = (x2-x1).^4+12.*x1.*x2-x1+x2-3;
% box on; mesh(x1,x2,f); %plotting the function 
contour(x1,x2,f, 20); hold on; 
for i = 1:length(X_init_array1)
    px(i) = X_init_array1{i}(1);
    py(i) = X_init_array1{i}(2);
end
hold on;
plot(px, py, '^-'); %plot sequence of points starting from starting point1

figure;
contour(x1,x2,f, 20); hold on; 
for i = 1:length(X_init_array2)
    px(i) = X_init_array2{i}(1);
    py(i) = X_init_array2{i}(2);
end
hold on;
plot(px, py, '^-'); %plot sequence of points starting from starting point1
%Steepest gradient
syms x1 x2;
f = (x2-x1).^4+12.*x1.*x2-x1+x2-3;
X_init_array1 = steepestGrad(f, X_init_array1);
X_init_array2 = steepestGrad(f, X_init_array2);