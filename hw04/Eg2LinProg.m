clc; close all;
f = [0 -10 0 -6 -20];
A = [];
b = [];
Aeq = [1 -1 -1 0 0; 0 0 1 -1 -1];
beq = [0;0];
lb = [0 0 0 0 0];
ub = [4 3 3 2 2];
[x,fval] = linprog(f, A, b, Aeq, beq, lb, ub);
for i = 1:length(x)
    fprintf('I%d = %d\n', i, x(i));
end
fprintf('The maximum value of function is %d\n', -fval);