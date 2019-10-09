function [] = steepestGrad(func,X_Init)
%UNTITLED6 Summary ofunc this funcunction goes here
%   Detailed explanation goes here
    syms x1 x2;
    epsilon = .001;
    dx1=diff(func,x1);
    grad{1}(1) = double(subs(dx1,{x1,x2},X_Init{1}));
    dx2=diff(func,x2);
    grad{1}(2) = double(subs(dx2,{x1,x2},X_Init{1}));
    X_Init{2}= X_Init{1}-alpha*grad{1};
    k = 2;    
    while(norm((X_Init{k-1}-X_Init{k}),2)>epsilon)
        dx1=diff(func,x1);
        grad{k}(1) = double(subs(dx1,{x1,x2},X_Init{k}));
        dx2=diff(func,x2);
        grad{k}(2) = double(subs(dx2,{x1,x2},X_Init{k}));
        X_Init{k+1}= X_Init{k}-alpha*grad{k};
        k = k+1;
    end
    x = linspace(-1,1,50);
    y = x;
    [x1,x2] = meshgrid(x,y);
    func = (x2-x1).^4+12.*x1.*x2-x1+x2-3;
    box on; mesh(x1,x2,func);
    funcigure;
    contour(x1,x2,func, 20); hold on; 
    for i = 1:length(X_Init)
        px(i) = X_Init{i}(1);
        py(i) = X_Init{i}(2);
    end
    hold on;
    plot(px, py, 'x-');
end