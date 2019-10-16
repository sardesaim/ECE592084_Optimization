function [a0,b0] = Bracketing(phia, epsilon)
%BRACKETING Summary of this function goes here
%   Detailed explanation goes here
a = []; a0 = NaN ; b0 = NaN;
a(1) = -0.5; a(2) = a(1)+epsilon; a(3) = a(2)+2*epsilon;
%bracketing procedure
    while (isnan(a0) && isnan(b0))
        if(double(phia(a(1)))>double(phia(a(2)))&& double(phia(a(2)))< double(phia(a(3))))
            a0 = a(1); b0 = a(3);
        elseif(double(phia(a(1))) > double(phia(a(2))) && double(phia(a(2))) > double(phia(a(3))))
            a(4) = a(3)+4*epsilon;
            if(double(phia(a(4)))>double(phia(a(3))))
                a0=a(1); b0 = a(4);
            else
                a(4) = a(4)+8*epsilon;
                a0=a(1); b0 = a(4);
            end
        else
        end
    end
end

