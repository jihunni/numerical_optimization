function [phenotype] = convert_to_solution_space(legnth_of_genotype, population)
%convert_to_solution_space Summary of this function goes here
% The converter of a vector from genotype space to phenotype space. For
% simplification, this program implementation assumes that search region is
% on [0, 1]
% (solution space) for genetic algorithm.
%   Detailed explanation goes here
%   Input : 
%            legnth_of_genotype : (integer) the dimension of genotype
%            termination condition : (integer) the number of generation as
%            a termination condition for optimization algorithm.
%            population_: (float vector) the population that will be
%            converted from genotype space to phenotype space.
%            
%   output : 
%            phenotype : (float vector) the converted vector in phenotype space.
converter = ones(1, legnth_of_genotype); % uniformly distributed on the search space
for k = 1:legnth_of_genotype
    converter(k) = converter(k) * (2^(k-1)) / (2^legnth_of_genotype);
end

phenotype = sum(population * transpose(converter) , 2);
%evaluation = model(convert_to_solution_space);
end

