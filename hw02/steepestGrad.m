function [X_init_array1] = steepestGrad(f, X_init_array1)
%GRADDESC Summary of this function goes here
%   Detailed explanation goes here
    syms x1 x2;
    epsilon = .001;
%     X_init_array1{1} = [-0.9, -0.5];
    %     Initialize first gradients and 2nd element
    grad{1} = double(subs(gradient(f(x1, x2), [x1 x2]), {x1,x2}, X_init_array1(1)))';
    syms alph;
    searchDir = -gradient(f(x1, x2), [x1 x2]);
    fx_ag = f(x1+alph*searchDir(1), x2+alph*searchDir(2));
    % fx_ag = subs(f(x1+alph*searchDir(1), x2+alph*searchDir(2)));
    f_alph = subs(fx_ag, [x1, x2], {X_init_array1(1)});
    f_dash_alph = diff(f_alph, alph);
    f_ddash_alph = diff(f_dash_alph, alph);
    alpha{1} = 0.02;
    alpha{2} = alpha{1}-double(subs((f_dash_alph/f_ddash_alph), alph, alpha{1}));
    X_init_array1{2}= X_init_array1{1}-alpha{1}*grad{1};
    k = 2;
    %     Iterate till the norm of consecutive points becomes less than epsilon
    while(k<500&& norm((X_init_array1{k}-X_init_array1{k-1}),2)>epsilon)
        syms alph;
        grad{k} = (double(subs(gradient(f(x1, x2), [x1 x2]), {x1,x2}, X_init_array1(k))))';
        %     i = 2;
        %     alpha_new{1} = alpha{k-1};
        %     alpha_new{2} = 0;
        %     while(abs(alpha_new{end}-alpha_new{end-1})>0.001)
        f_alph = subs(fx_ag, [x1, x2], {X_init_array1(k)});
        f_dash_alph = diff(f_alph, alph);
        f_ddash_alph = diff(f_dash_alph, alph);
        alpha{k} = alpha{k-1}- double(subs((f_dash_alph/f_ddash_alph), alph, alpha{k-1}));
        %         alpha{k+1} = alpha{k}- double(subs((f_dash_alph/f_ddash_alph), alph, alpha{k}));
        %         i = i+1;
        %     end
        %     alpha{k} = alpha_new{i-1};
        X_init_array1{k+1}= X_init_array1{k}-alpha{k}*grad{k};
        k = k+1;
    end
end

