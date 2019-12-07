function [num] = bit2num(bit, range)
%BIT2NUM convert the binary encoded bits to actual numbers 
%   Detailed explanation goes here
    integer = polyval(bit, 2);
    num = range(1) + integer*(range(2)-range(1))/(2^length(bit)-1);
end


