function [fitness] = evalPopu(popu,bit_n,range,obj_func)
%EVALPOPU Evaluates fitness function
%   Split the binary encoded string into two parts and then evaluates the
%   fitness at that points. 
    [popu_s, string_length] = size(popu); 
    fitness = zeros(popu_s,1);
    for i = 1:popu_s
        num1 = bit2num(popu(i,1:bit_n), range);
        num2 = bit2num(popu(i,bit_n+1:string_length), range);
        fitness(i) = obj_func(num1, num2);
    end
end

