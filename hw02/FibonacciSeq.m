function [alpha] = FibonacciSeq(a0, b0, phia)
%FIBONACCISEQ Summary of this function goes here
%   Detailed explanation goes here
%     rho = 
    epsi = 0.075;
    N = ceil(1+2*0.01/(0.01/(b0-a0)));
    n=1;
    while (fibonacci(n)<N)
        fiboNum(n) = fibonacci(n+1);
        n=n+1;
    end
    a = a0;
    b = b0;
    rho(1) = 1-fiboNum(end)/fiboNum(end-1);
    s = a+rho(1)*(b-a);
    t = a+(1-rho(1))*(b-a);
    f1=double(phia(s));
    f2=double(phia(t));
    dat = {'Iteration','rhok','ak','bk','f(ak)','f(bk)','new int'};
    for n = 1:N
        rho(n) = 1-fiboNum(end-n)/fiboNum(end-n-1);
        if double(phia(s))< double(phia(t))
           b = t;
           t = s;
           s = a+rho(n)*(b-a);
           f2 = f1;
           f1 = double(phia(s));
        else 
           a = s;
           s = t;
           s = a+(1-rho(n))*(b-a);
           f1 = f2;
           f2 = double(phia(t));
        end
%         dat{n,:} = {n,rho,s,t,f1,f2,[s,t]};
    end
    alpha = (s+t)/2;
end
end

