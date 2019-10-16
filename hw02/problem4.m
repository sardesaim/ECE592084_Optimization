% problem4;
clear all; clc; close all;
% Use newtons method to minimize the given function
syms x1 x2;
%define the function in terms of symbolic variables x1 and x2
f = (x2-x1).^4+12.*x1.*x2-x1+x2-3;
X_init_array1{1} = [.4,-.45];
X_init_array1{2} = [0,0]; %initialize second point in the sequence
% alpha = 0.02;
epsilon = .001; %set epsilon to some small value
%calculate the gradient and hessian at the initial point x0
grad{1} = (double(subs(gradient(f, [x1 x2]), {x1,x2}, X_init_array1(1))));
hess{1} = double(subs(hessian(f, [x1 x2]), {x1,x2}, X_init_array1(1)));
%The conditions for the positive definiteness can be removed to check
%the updated sequence we get through the newtons method. But as this
%example has no pd matrix the sequence obtained by using the Hessian as the
%step size does not converge to the minimizer.
if (eig(hess{1})>0) %check for the kth Hessian to be pd
    %     update point 2 with hess{1} as step size
    X_init_array1{2}= X_init_array1{1}-(pinv(hess{1})*grad{1})';
    k = 2;
    %loop to update the points while calculating the hessian at every step
    %to set the step size
    while(norm((X_init_array1{k-1}-X_init_array1{k}),2)~=0)
        grad{k} = (double(subs(gradient(f, [x1 x2]), {x1,x2}, X_init_array1(k))));
        hess{k} = double(subs(hessian(f, [x1 x2]), {x1,x2}, X_init_array1(k)));
        if eig(hess{k})>0 %check if hess{k} is pd
            %             update kth point with hess{k} as step size
            X_init_array1{k+1}= X_init_array1{k}-(pinv(hess{k})*grad{k})';
            k = k+1;
        else
            disp('Hessian of the matrix is not pd. Newtons method would not converge');
        end
    end
else
    disp('Hessian of the matrix is not pd. Newtons method would not converge');
end

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
plot(px, py, 'x-'); %plot sequence of points starting from starting point1
