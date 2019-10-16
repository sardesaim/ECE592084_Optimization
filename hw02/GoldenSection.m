function [alpha,dat] = GoldenSection(a0, b0, phia)
%GOLDENSECTION Summary of this function goes here
%   Detailed explanation goes here
%golden section
    epsi = 0.075;
    N = ceil(log(0.01/(b0-a0))/log(0.61803));
    rho = 0.382;
    a = a0;
    b = b0;
    s = a+rho*(b-a);
    t = a+(1-rho)*(b-a);
    f1=double(phia(s));
    f2=double(phia(t));
    dat = {'Iteration','rhok','ak','bk','f(ak)','f(bk)','new int'};
    for n = 1:N
        if double(phia(s))< double(phia(t))
           b = t;
           t = s;
           s = a+rho*(b-a);
           f2 = f1;
           f1 = double(phia(s));
        else 
           a = s;
           s = t;
           s = a+(1-rho)*(b-a);
           f1 = f2;
           f2 = double(phia(t));
        end
%         dat{n,:} = {n,rho,s,t,f1,f2,[s,t]};
    end
    alpha = (s+t)/2;
end

