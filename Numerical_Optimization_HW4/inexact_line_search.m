function [alpha] = inexact_line_search(objective_function,x, y, P, objective_function_gradient)
%INEXACT_LINE_SEARCH Summary of this function goes here
%Inesact line search is an algorithm that calculates a proper alpha
%(step length) for gradient-based method in heuristic way.
%   Detailed explanation goes here
%   Input : 
%            objective_function
%            x
%            y
%            P : 2 * 1 vector 
%            objective_funobejctivection_gradient : 2 * 1 vector
%   output : 
%            alpha
    
    % to measure the performance
    %tic
    iteration = 0;
    
    % initialize
    alpha = random('Uniform', 80,100) ;
    %contraction_factor = random('Uniform', 0.1,1) ;
    contraction_factor = 0.9 ;
    wolfe_constant = random('Uniform', 0.1,1) ;
    
    while( ...
            objective_function(x + alpha*P(1), y + alpha*P(2)) > ...
            alpha * wolfe_constant * (objective_function_gradient * transpose(P)) + objective_function(x, y) & ...
            alpha > 0.1 ...
        )
        alpha = alpha * contraction_factor ;
        iteration = iteration+1;
    end
    %fprintf('alpha: %i\n', alpha)
    %fprintf('the number of iteration: %i\n', iteration)
    %toc
end

