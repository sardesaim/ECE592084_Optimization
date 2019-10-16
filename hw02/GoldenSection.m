function [s,t, dat] = GoldenSection(a0, b0, fcn)
%GOLDENSECTION Summary of this function goes here
%   Detailed explanation goes here
%golden section
    epsi = 0.075;
    N = floor(log(0.01/norm(b0-a0))/log(0.6180));
    rho = 0.382;
    a = a0;
    b = b0;
    s = a+rho*(b-a);
    t = a+(1-rho)*(b-a);
    f1=double(fcn(s(1), s(2)));
    f2=double(fcn(t(1), t(2)));
    dat = {'Iteration','rhok','ak','bk','f(ak)','f(bk)','new int'};
    for n = 1:N
        if double(fcn(s(1), s(2)))< double(fcn(t(1), t(2)))
           b = t;
           t = s;
           s = a+rho*(b-a);
           f2 = f1;
           f1 = double(fcn(s(1), s(2)));
        else 
           a = s;
           s = t;
           s = a+(1-rho)*(b-a);
           f1 = f2;
           f2 = double(fcn(t(1), t(2)));
        end
%         dat{n,:} = {n,rho,s,t,f1,f2,[s,t]};
        dat{n+1,1} = n;
        dat{n+1,2} = rho;
        dat{n+1,3} = s;
        dat{n+1,4} = t;
        dat{n+1,5} = f1;
        dat{n+1,6} = f2;
        dat{n+1,7} = [s,t];
    end
end

