% problem4syms x1 x2;
syms x1 x2;
f = (x2-x1)^4+12*x1*x2-x1+x2-3;
X_init_array1{1} = [0.55, 0.7];
X_init_array1{2} = [0,0];

alpha = 0.02;
epsilon = .001;

dx1=diff(f,x1);
g{1}(1) = double(subs(dx1,{x1,x2},X_init_array1{1}));
dx2=diff(f,x2);
g{1}(2) = double(subs(dx2,{x1,x2},X_init_array1{1}));
X_init_array1{2}= X_init_array1{1}-alpha*g{1};