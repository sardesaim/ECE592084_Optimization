function [X_init_array1] = RankOne(f, X_init_array1)
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
    H{1} = [2 1; 1 2];
    searchDir = -gradient(f(x1, x2), [x1 x2]); %search dir in direction of negative grad
    k = 1;
    dk = @(x1, x2,k) (-H{k}*searchDir);
    fx_ag = f(x1+alph*dk(x1,x2,k), x2+alph(x1,x2)*dk(x1, x2,k)); %f(x-alpha*g)
    f_alph = subs(fx_ag, [x1, x2], {X_init_array1(1)}); %f(x1-alph*g1)
    [s,t,~] = GoldenSection(0,0.5,f_alph);
    alpha{1} = (s+t)/2;
    %find 2nd point based on initial alpha
    X_init_array1{2}= X_init_array1{1}-alpha{1}*dk(X_init_array1{1}(1), X_init_array1{1}(2));
    deltaX{1} = alpha{1}*dk;
    deltaG{1} = double(subs(gradient(f(x1, x2), [x1 x2]), {x1,x2}, X_init_array1(2)))'...
        -double(subs(gradient(f(x1, x2), [x1 x2]), {x1,x2}, X_init_array1(1)))';
    H{2} = H{1} + ((deltaX{1}-H{1}*deltaG{1})*(deltaX{1}-H{1}*deltaG{1})')/(...
        deltaG{1}'*(deltaX{1}-H{1}*deltaG{1}));
    k = 2;
    %     Iterate till the norm of consecutive points becomes less than epsilon
    while(k<100 && norm((X_init_array1{k}-X_init_array1{k-1}),2)>epsilon)
        syms alph;
        dk = -H{k}*searchDir;
        grad{k} = (double(subs(gradient(f(x1, x2), [x1 x2]), {x1,x2}, X_init_array1(k))))';
        f_alph = subs(fx_ag, [x1, x2], {X_init_array1(k)});
        [s,t,~] = GoldenSection(0,0.5, f_alph);
        alpha{k} = (s+t)/2;
        X_init_array1{k+1}= X_init_array1{k}-alpha{k}*dk{k}; %update point using the new found alpha
        deltaX{k} = X_init_array1{k}-X_init_array1{k};
    deltaG{1} = double(subs(gradient(f(x1, x2), [x1 x2]), {x1,x2}, X_init_array1(2)))'...
        -double(subs(gradient(f(x1, x2), [x1 x2]), {x1,x2}, X_init_array1(1)))';
    H{2} = H{1} + ((deltaX{1}-H{1}*deltaG{1})*(deltaX{1}-H{1}*deltaG{1})')/(...
        deltaG{1}'*(deltaX{1}-H{1}*deltaG{1}));
        k = k+1;
    end
end