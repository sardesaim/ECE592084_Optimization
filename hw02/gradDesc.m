function [X_init_array1] = gradDesc(f, X_init_array1)
%GRADDESC Summary of this function goes here
%   Detailed explanation goes here
    syms x1 x2;
    grad = {};
    alpha = 0.02;
    epsilon = .001;
%     Initialize first gradients and 2nd element
    grad{1} = (double(subs(gradient(f, [x1 x2]), {x1,x2}, X_init_array1(1))))';
    X_init_array1{2}= X_init_array1{1}-alpha*grad{1};
    k = 2;    
%     Iterate till the norm of consecutive points becomes less than epsilon
    while(norm((X_init_array1{k-1}-X_init_array1{k}),2)>epsilon)
        grad{k} = (double(subs(gradient(f, [x1 x2]), {x1,x2}, X_init_array1(k))))';
        X_init_array1{k+1}= X_init_array1{k}-alpha*grad{k};
        k = k+1;
    end
end

