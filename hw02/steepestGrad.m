function [X_init_array1] = steepestGrad(f, X_init_array1)
%GRADDESC Summary of this function goes here
%   Detailed explanation goes here
    syms x1 x2 alph;
    grad = {};
    epsilon = .001;
%     Initialize first gradients and 2nd element
    grad{1} = subs(gradient(f, [x1 x2]), {x1,x2}, X_init_array1(1))';
    searchDir = -gradient(f, [x1 x2]);
    fx_ag = subs(f, [x1, x2], {x1+alph*searchDir(1), x2+alph*searchDir(2)});
    f_alph{1} = subs(fx_ag, [x1, x2], {X_init_array1(1)});
    f_dash_alph{1} = diff(f_alph{1}, alph);
    f_ddash_alph{1} = diff(diff(f_alph{1}, alph));
    alpha{1} = 0.05;
    alpha{2} = alpha{1} - double(subs(f_dash_alph{1}/f_ddash_alph{1}, alph, alpha{1}));
    i = 2
%     new_f = subs(f-alph*gradient(f,[x1 x2]),{x1 x2}, [X_init_array1(1)]);
%     new_fn = subs(f,new_f);
%     sub
%     alpha{1} = argmin(double(subs(f-alph*gradient(f,[x1 x2]),{x1 x2}, [X_init_array1(1)])));
%     X_init_array1{2}= X_init_array1{1}-alpha{1}*grad{1};
    X_init_array1{2}= X_init_array1{1}-alpha{1}*grad{1};
    k = 2;    
%     Iterate till the norm of consecutive points becomes less than epsilon
    while(norm((X_init_array1{k-1}-X_init_array1{k}),2)>epsilon)
        syms alph;
        g = subs(f, [x1 x2], [X_init_array1(k)+alph*searchDir]);
        grad{k} = (double(subs(gradient(f, [x1 x2]), {x1,x2}, X_init_array1(k))))';
        while(norm(alpha{i-1}-alpha(i)>epsilon))
            f_alph{i} = subs(fx_ag, [x1, x2], {X_init_array1(1)});
            f_dash_alph{i} = diff(f_alph, alph);
            f_ddash_alph{i} = diff(diff(f_alph, alph)); 
            alpha{i+1} = alpha{i} - double(subs(f_dash_alph{1}/f_ddash_alph{1}, alph, alpha{i}));
        end
        X_init_array1{k+1}= X_init_array1{k}-alpha{i}*grad{k};
        k = k+1;
    end
end

