function [s,t, dat] = FibonacciSeq(a0, b0, fcn)
%FIBONACCISEQ Summary of this function goes here
%   Detailed explanation goes here
%     rho = 
    N = ceil(1+(2*0.1)/(0.01/norm(b0-a0)));
    n=1;
    while (fibonacci(n)<N)
        fiboNum(n) = fibonacci(n+1);
        n=n+1;
    end
    for i = 1:length(fiboNum)-1
        rho(i) = 1- (fiboNum(length(fiboNum)-i)/fiboNum(length(fiboNum)-i+1));
    end
    a = a0;
    b = b0;
    s = a+rho(1)*(b-a);
    t = a+(1-rho(1))*(b-a);
    f1= fcn(s(1), s(2));
    f2= fcn(t(1), t(2));
    dat = {'Iteration','rhok','ak','bk','f(ak)','f(bk)','new int'};
    for n = 1:length(fiboNum)-1
        if fcn(s(1), s(2))< fcn(t(1), t(2))
           b = t;
           t = s;
           s = a+rho(n)*(b-a);
           f2 = f1;
           f1 = fcn(s(1), s(2));
        else 
           a = s;
           s = t;
           t = a+(1-rho(n))*(b-a);
           f1 = f2;
           f2 = fcn(t(1), t(2));
        end
%         dat{n,:} = {n,rho,s,t,f1,f2,[s,t]};
        dat{n+1,1} = n;
        dat{n+1,2} = rho(n);
        dat{n+1,3} = mat2str(s);
        dat{n+1,4} = mat2str(t);
        dat{n+1,5} = f1;
        dat{n+1,6} = f2;
        dat{n+1,7} = mat2str([s;t]);
    end
end