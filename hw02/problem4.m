% problem4syms x1 x2;
syms x1 x2;
f = (x2-x1).^4+12.*x1.*x2-x1+x2-3;
X_init_array1{1} = [0.55,0.7];
X_init_array1{2} = [0,0];

alpha = 0.02;

epsilon = .001;
grad{1} = (double(subs(gradient(f, [x1 x2]), {x1,x2}, X_init_array1(1))));
hess{1} = double(subs(hessian(f, [x1 x2]), {x1,x2}, X_init_array1(1)));
X_init_array1{2}= X_init_array1{1}-(pinv(hess{1})*grad{1})';

 k = 2;    
 while(norm((X_init_array1{k-1}-X_init_array1{k}),2)~=0)
    grad{k} = (double(subs(gradient(f, [x1 x2]), {x1,x2}, X_init_array1(k))));
    hess{k} = double(subs(hessian(f, [x1 x2]), {x1,x2}, X_init_array1(k)));
    X_init_array1{k+1}= X_init_array1{k}-(pinv(hess{k})*grad{k})';
    k = k+1;
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
