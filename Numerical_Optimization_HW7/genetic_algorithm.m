function [X,evaluation_score] = genetic_algorithm(model, legnth_of_genotype, termination_iteration, population_size, mutation_probability, crossover_probability)
%GENETIC_ALGORITHM Summary of this function goes here
% Univariable global stochastic optimization algorithm. For simplification,
% this program search only on [0, 1]
%   Detailed explanation goes here
%   Input : 
%            model : (univariate objective function) to minimize
%            legnth_of_genotype : (integer) the dimension of genotype
%            termination condition : (integer) the number of generation as
%            a termination condition for optimization algorithm.
%            population_size : (integer) the number of population for
%            genetic algorithmn.
%            mutation_probability : (float, [0, 1]) the mutation
%            probability of parent to reproduce children generation (e.g.
%            0.1)
%            corssover_probability : (float, [0,1]) the one-point crossover
%            probability of parent to reproduce children generation. (e.g.
%            80~90 %)
%            
%   output : 
%            x : (float) global minimuma on [0, 1]
%            evoluation_score : (float) global minimum function value
%   function :
%            convert_to_solution_space : the function that converts a vecotr from
%            genetic space to solution sapce
% to measure the performance
tic

% initialize a population of chromosomes
num_iteration = 0;
population = randi([0 1], population_size, legnth_of_genotype);
    % each row corresponds to one sample.

% evaluation the fitness of each individual : skip
% to measure the performance
    tic

% termination criterion
while num_iteration < termination_iteration
    % select parents
    evaluation_score = (-1) * model(convert_to_solution_space(legnth_of_genotype, population)) ;
    selection_probability = (evaluation_score-min(evaluation_score)) ./ (max(evaluation_score) - min(evaluation_score)) ;
    selection = rand(population_size, 1) < selection_probability ;
    parent = population(find(selection), :);
    clear selection

    % reproduce new individuals
    % (1) mutation
    boolian_mutation = find(rand(size(parent, 1), 1) < mutation_probability) ;
    for j = boolian_mutation
        if parent(j,randi(legnth_of_genotype)) == 0
            parent(j,randi(legnth_of_genotype)) = 1 ;
        else
            parent(j,randi(legnth_of_genotype)) = 0 ;
        end
    end
    clear boolian_mutation
    % (2) one point cross-over
    boolian_crossover = find(rand(size(parent, 1), 1) < crossover_probability) ;
    % error in j
    for j = transpose(boolian_crossover)
        %fprintf("%d", j)
        another_parent = randi(size(parent,1));
        corss_over_point = randi(legnth_of_genotype);
        tmp = parent(j, 1:corss_over_point) ; %temporary 
        parent(j, 1:corss_over_point) = parent(another_parent, 1:corss_over_point) ;
        parent(another_parent, 1:corss_over_point) = tmp ;
    end
    children = parent;
    clear parent;
    clear boolian_crossover ;
    clear tmp ;
    clear corss_over_point ;
    clear another_parent ;

    % evaluation the fitness of each new individuals
    evaluation_score = model(convert_to_solution_space(legnth_of_genotype, children)) ;
    %fprintf("evaluation score : %f\n", evaluation_score)
    
    % replace into new individuals (randomly)
    selection = randi(size(population, 1), size(population, 1) - size(children, 1),1);
    population = [population(selection,:); children]; % replace into new population

    num_iteration = num_iteration + 1 ;
end

X = convert_to_solution_space(legnth_of_genotype, population) ;
[evaluation_score, index] = min(model(X)) ;
X = X(index);

fprintf("optimal X : %f\n", X)
fprintf("optimal value : %f\n", evaluation_score)
toc
end

