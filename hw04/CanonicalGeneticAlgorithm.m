clc; close all;
%% Init
generation_n = 100; 
popuSize = 100; 
xover_rate = 0.75; 
mutate_rate = 0.01;
%standard form of canonical GA does maximization hence '-' sign to the
%objective function. 
obj_func = @(x1, x2) -((x2-x1).^4 + 12.*x1.*x2 - x1 + x2 - 3);
var_n = 2;
range = [-1 1];
init_popu = rand(popuSize,2)*2.0-1; %initialize random population
%convert the random initialization into encoded binary string
r = 10^-4;  %resolution 
%calculate the number of bits required to represent the float in a binary string
li = ceil(log2((max(init_popu)-min(init_popu)+r)/r));   
% li = ceil(log2((1-(-1)+r)/r));
bit_n = sum(li);    %calculate total number of bits 
n = ceil(bit_n/8);   % number bits for integer part of your number
m = bit_n-n;         % number bits for fraction part of your number
% binary number conversion for the intitial population 
popu = [fix(rem(init_popu(:,1).*pow2(-(n-1):m),2)),...
    fix(rem(init_popu(:,2).*pow2(-(n-1):m),2))];
% popu = rand(popuSize, bit_n * var_n)>0.5;
upper = zeros(generation_n, 1); 
average = zeros(generation_n, 1);
lower = zeros(generation_n, 1);
%% Run for n generations. 
for i = 1:generation_n
    func_value = evalPopu(popu, bit_n, range, obj_func);
%     fill obj func matrices
    upper(i) = max(func_value); 
    average(i) = mean(func_value); 
    lower(i) = min(func_value);
    %calculate the next population 
    popu = nextPopu(popu, func_value, xover_rate, mutate_rate);
end
% Convert Final population binary encoded to float numbers
final_popu = zeros(popuSize,2);
for i = 1:popuSize
    final_popu(i,1) = bit2num(popu(i,1:bit_n), range);
    final_popu(i,2) = bit2num(popu(i,bit_n+1:end), range);
end
[~,index] = max(obj_func(final_popu(:,1), final_popu(:,2)));
%% plotting 
epochs = 1:popuSize;
X = linspace(-1,1,1000);
Y = X;
[Xpt, Ypt] = meshgrid(X,Y);
Z = obj_func(Xpt,Ypt);
figure;
contour(Xpt,Ypt, Z, 30); hold on;
plot(final_popu(index,1), final_popu(index,2), 'x');
figure;
plot(epochs, upper, 'r', epochs, average, 'b', epochs, lower, 'g');
legend('best','average','worst');