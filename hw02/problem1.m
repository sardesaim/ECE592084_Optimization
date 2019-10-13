clc; clear all; close all;

syms x1 x2 alph;
f = @(x1, x2) x1.^2+x2.^2+x1.*x2;
X_init{1} = [0.8,0.25];
gx = @(f) gradient(f(x1, x2), [x1 x2])';
a = 0.3;
X_init{2} = X_init{1} - a*double(subs(gx(f)),[x1 x2], X_init(1)));
phia = f(x1-alph*gx(f(x1, x2)), x2-alph*gx(f(x1, x2)));
