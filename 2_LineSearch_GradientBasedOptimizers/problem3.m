clc; clear all; close all;

%Gradient descent 
syms x1 x2;
f = @(x1, x2) (x2-x1).^4+12.*x1.*x2-x1+x2-3; %declare the function in terms of x1, x2
%Initializations for the sequence of points
X_init_array1{1} = [0.55, 0.7]; %starting point 1
X_init_array1{2} = [0,0];
X_init_array2{1} = [-0.9 -0.5]; %starting point 2
X_init_array2{2} = [0,0];
%Gradient descent from the two starting points 
X_init_array1 = gradDesc(f, X_init_array1);
X_init_array2 = gradDesc(f, X_init_array2);

%Steepest gradient method
X_init_array3{1} = [0.55, 0.7]; %starting point 1
X_init_array3{2} = [0,0];
X_init_array3 = steepestGrad(f, X_init_array3);

X_init_array4{1} = [-0.9,-0.5]; %starting point 1
X_init_array4{2} = [0,0];
X_init_array4 = steepestGrad(f, X_init_array4);

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
    px1(i) = X_init_array2{i}(1);
    py1(i) = X_init_array2{i}(2);
end
hold on;
plot(px1, py1, 'x-'); %plot sequence of points starting from starting point2

figure;
contour(x1,x2,f, 20); hold on; 
for i = 1:length(X_init_array3)
    px2(i) = X_init_array3{i}(1);
    py2(i) = X_init_array3{i}(2);
end
hold on;
plot(px2, py2, 'x-'); %plot sequence of points starting from starting point1 using sd

figure;
contour(x1,x2,f, 20); hold on; 
for i = 1:length(X_init_array4)
    px3(i) = X_init_array4{i}(1);
    py3(i) = X_init_array4{i}(2);
end
hold on;
plot(px3, py3, 'x-'); %plot sequence of points starting from starting point1 using sd