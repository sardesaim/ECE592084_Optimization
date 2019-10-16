function [a0,b0] = Bracketing(f, x_init, epsilon)
%BRACKETING Summary of this function goes here
%   Detailed explanation goes here
   x0 = x_init; x1 = x0+epsilon; x2 = x1+2*epsilon;
    n = 1;
    %bracketing procedure
    if((f(x0(1), x0(2))>f(x1(1), x1(2))) && (f(x1(1), x1(2))<f(x2(1), x2(2))))
            a0 = x0; 
            b0 = x2;
    else
        while not((f(x0(1), x0(2))>f(x1(1), x1(2))) && (f(x1(1), x1(2))<f(x2(1), x2(2))))
            if((f(x0(1), x0(2))>f(x1(1), x1(2))) && (f(x1(1), x1(2))>f(x2(1), x2(2))))
                x0 = x1; 
                x1 = x2; 
                x2 = x2+(2^(n+1))*epsilon;
            elseif((f(x0(1), x0(2))<f(x2(1), x2(2))) && (f(x1(1), x2(2))<f(x2(1), x2(2))))
                x2 = x1;
                x1 = x0;       
                x0 = x0-(2^(n+1))*epsilon;
            end
            n = n+1;
        end
        a0 = x0; 
        b0 = x2;
    end
end