function [s,t, dat] = GoldenSection(a0, b0, fcn)
%GOLDENSECTION Summary of this function goes here
%   Detailed explanation goes here
%golden section
    epsi = 0.075;
    N = ceil(log(0.01/norm(b0-a0))/log(0.6180));
    rho = 0.382;
    a = a0;
    b = b0;
    s = a+rho*(b-a);
    t = a+(1-rho)*(b-a);
    f1= fcn(s(1), s(2));
    f2= fcn(t(1), t(2));
    dat = {'Iteration','rhok','ak','bk','f(ak)','f(bk)','new int'};
    for n = 1:N
        if fcn(s(1), s(2))< fcn(t(1), t(2))
           b = t;
           t = s;
           s = a+rho*(b-a);
           f2 = f1;
           f1 = fcn(s(1), s(2));
           dat{n+1,5} = f1;
           dat{n+1,6} = f2;    
        elseif fcn(s(1), s(2))> fcn(t(1), t(2))
           a = s;
           s = t;
           t = a+(1-rho)*(b-a);
           f1 = f2;
           f2 = fcn(t(1), t(2));
           dat{n+1,5} = f1;
           dat{n+1,6} = f2;
        else
            break;
        end
%         dat{n,:} = {n,rho,s,t,f1,f2,[s,t]};
        dat{n+1,1} = n;
        dat{n+1,2} = rho;
        dat{n+1,3} = mat2str(s);
        dat{n+1,4} = mat2str(t);
        dat{n+1,7} = mat2str([s;t]);
    end
end

