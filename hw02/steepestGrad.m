function [X_init_array1] = steepestGrad(f, X_init_array1)
%Steepest Gradient - Used to calculate the minima using the steepest
%gradient method
%   The function takes f - The objective function and the Initial point
%   from the seqeuence as the input and gives out an array of sequence of
%   points to reach the minima as output. At each step Line search is used
%   to find the optimum step size
    syms x1 x2;
    epsilon = .001;
%     Initialize first gradients and 2nd element
    grad{1} = double(subs(gradient(f(x1, x2), [x1 x2]), {x1,x2}, X_init_array1(1)))';
    syms alph;
    searchDir = -gradient(f(x1, x2), [x1 x2]); %search dir in direction of negative grad
    fx_ag = f(x1+alph*searchDir(1), x2+alph*searchDir(2)); %f(x-alpha*g)
    f_alph = subs(fx_ag, [x1, x2], {X_init_array1(1)}); %f(x1-alph*g1)
    f_dash_alph = diff(f_alph, alph); %f'(x1-alph*g1)
    f_ddash_alph = diff(f_dash_alph, alph); %f''(x1-alph*g1)
    %initialize alpha values to arbitrary levels
    alpha{1} = 0.02;
    alpha{2} = alpha{1}-double(subs((f_dash_alph/f_ddash_alph), alph, alpha{1}));
    %find 2nd point based on initial alpha
    X_init_array1{2}= X_init_array1{1}-alpha{1}*grad{1};
    k = 2;
    %     Iterate till the norm of consecutive points becomes less than epsilon
    while(k<100&& norm((X_init_array1{k}-X_init_array1{k-1}),2)>epsilon)
        syms alph;
        %find the kth gradient, f(xk-alph*gk), f'(xk-alph*gk),
        %f''(xk-alph*gk)
        grad{k} = (double(subs(gradient(f(x1, x2), [x1 x2]), {x1,x2}, X_init_array1(k))))';
        f_alph = subs(fx_ag, [x1, x2], {X_init_array1(k)});
        f_dash_alph = diff(f_alph, alph);
        f_ddash_alph = diff(f_dash_alph, alph);    
        %initialize alpha for line search using NM
        alpha_new{1} = alpha{k-1};
        alpha_new{2} = alpha{k};
        i = 2;
%         while loop to evalulate optimum alpha to minimize f(x-alpha*g)
%         using Newtons method (Line search)
        while(norm(alpha_new{i}-alpha_new{i-1})>0.001)
%         this while loop can be removed and the statement below can be
%         used to get alpha after one iteration of NM. And the point
%         converges to minima using this approach. Not sure if it is the
%         correct way. 
%             alpha{k} = alpha{k-1}- double(subs((f_dash_alph/f_ddash_alph), alph, alpha{k-1}));
            alpha_new{i+1} = alpha_new{i}- double(subs((f_dash_alph/f_ddash_alph), alph, alpha_new{i}));
            i = i+1;
        end
        alpha{k+1} = alpha_new{i}; %update alpha with the one found using NM
        X_init_array1{k+1}= X_init_array1{k}-alpha{k+1}*grad{k}; %update point using the new found alpha
        k = k+1;
    end
end
