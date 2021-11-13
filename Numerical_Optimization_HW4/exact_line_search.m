function [alpha] = exact_line_search(objective_function,x, y, P, objective_function_gradient)
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
    wolfe_constant_1 = random('Uniform', 0.1,1) ;
    wolfe_constant_2 = random('Uniform', 0.1, wolfe_constant_1);
    grid_interval = 0.001 ;
    objective_function_gradient_nextPoint = ...
        [(objective_function(x + alpha*P(1)+grid_interval, y + alpha*P(2)) - objective_function(x + alpha*P(1), y + alpha*P(2)))/grid_interval; ...
        (objective_function(x + alpha*P(1), y + alpha*P(2)+grid_interval) - objective_function(x + alpha*P(1), y + alpha*P(2)))/grid_interval];

    while( ...
            objective_function(x + alpha*P(1), y + alpha*P(2)) > ...
            wolfe_constant_1 * alpha * transpose(objective_function_gradient) * P + objective_function(x, y) | ...
            transpose(objective_function_gradient_nextPoint) * P < wolfe_constant_2 * transpose(objective_function_gradient) * P ...
        )
        alpha = alpha * contraction_factor ;
        iteration = iteration+1;
        objective_function_gradient_nextPoint = ...
        [(objective_function(x + alpha*P(1)+grid_interval, y + alpha*P(2)) - objective_function(x + alpha*P(1), y + alpha*P(2)))/grid_interval; ...
        (objective_function(x + alpha*P(1), y + alpha*P(2)+grid_interval) - objective_function(x + alpha*P(1), y + alpha*P(2)))/grid_interval];

        if (alpha < 0.001)
            fprintf("break by low alpha\n")
            break
        end
    end
    fprintf('alpha: %i\n', alpha)
    fprintf('the number of iteration: %i\n', iteration)
    %toc
end

