clc; clear all; close all;

% Rank One Correction
syms x1 x2;
f = @(x1, x2) (x2-x1).^4+12.*x1.*x2-x1+x2-3; %declare the function in terms of x1, x2
%Initializations for the sequence of points
X_init_array1{1} = [0.55, 0.7]; %starting point 1
X_init_array1{2} = [0,0];
X_init_array2{1} = [-0.9 -0.5]; %starting point 2
X_init_array2{2} = [0,0];
