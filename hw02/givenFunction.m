function [val] = givenFunction(x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    x1 = x{1}(1);
    x2 = x{1}(2);
    val = (x2-x1)^4+12*x1*x2-x1-x2-3;
end

