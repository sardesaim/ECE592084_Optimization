function [new_popu] = nextPopu(popu, fitness, x_over,mut_rate)
%NEXTPOPU Generate the next generation population 
%   Use the roulette wheel algorithm for selection and then does crossover
%   and mutation with the supplied probabilities to the function. 
    new_popu = popu;
    [popu_s, string_length] = size(popu); 
    %rescaling the fitness
    fitness = fitness-min(fitness);
    %find the best two and keep them 
    tmp_fitness = fitness; 
    [~, index1] = max(tmp_fitness); %find best 
    tmp_fitness(index1) = min(tmp_fitness); 
    [~, index2] = max(tmp_fitness); %find 2nd best 
    new_popu([1 2], :) = popu([index1 index2],:);
    %for roulette wheel algo 
    total = sum(fitness);
    if total == 0 
        fitness = ones(popu_s, 1)/popu_s; %sum of probs is 1
    end
    cumprob = cumsum(fitness)/total;
    %selection and crossover;
    for i = 2:popu_s/2
        tmp = find(cumprob - rand>0);
        parent1 = popu(tmp(1), :);
        tmp = find(cumprob - rand>0); 
        parent2 = popu(tmp(1), :);
        %do crossover 
        if rand<x_over
            xover_pt = ceil(rand*(string_length-1));
            new_popu(i*2-1,:) = [parent1(1:xover_pt) parent2(xover_pt+1:string_length)];
            new_popu(i*2,:) = [parent2(1:xover_pt) parent1(xover_pt+1:string_length)];
        end
    end
    %mutation 
    mask = rand(popu_s, string_length)<mut_rate;
    new_popu = xor(new_popu,mask);
    %restore elites 
    new_popu([1 2], :) = popu([index1 index2],:);
end