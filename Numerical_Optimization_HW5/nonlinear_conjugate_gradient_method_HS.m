function [x, y, x_all_iteration, y_all_iteration] = nonlinear_conjugate_gradient_method_HS(objective_function, x, y, termination_condition)
%nonlinear_conjugate_gradient_method_Hestenes-Stiefel(CG-HR) Summary of this function goes here
% Multivariate optimization algorithm ; Netwon's method is a gradient-based
% optimization method. This program is implemented for 2 dimensional
% objective function f(x,y)
%   Detailed explanation goes here
%   Input : 
%            objective_function : (function) to minimize
%            x : (float) initial point
%            y : (float) inital poin
%            termination condition : (float) termination condition for
%            optimization algorithm
%   output : 
%            x : (float) local minimuma
%            y : (float) local minimuma
%            x_all_iteration : (list) iterative history of x (for plot)
%            y_all_iteration : (list) iterative history of y (for plot)
 % to measure the performance
    tic
    iteration = 0;
    % initialization
    x_all_iteration=[x];
    y_all_iteration=[y];
    grid_size = 0.001;
    gradient = [ (objective_function(x+grid_size, y) - objective_function(x, y))/ grid_size ; (objective_function(x, y+grid_size) - objective_function(x, y))/ grid_size] ;
    P = -gradient;

    %while(iteration < 5 | ...
    %        sqrt((x-x_all_iteration(end-1))^2 + (y-y_all_iteration(end-1)^2)) > termination_condition ...
    %      )
    while(transpose(gradient) * gradient >= termination_condition)
        %alpha = -(transpose(gamma) * P) / (transpose(P) * linear_A * P) ;
        alpha =  exact_line_search_with_strong_Wolfe_condition(objective_function, x, y, P, gradient);
        x = x + alpha * P(1) ;
        y = y + alpha * P(2) ;
        gradient_old = gradient;
        gradient = [ (objective_function(x+grid_size, y) - objective_function(x, y))/ grid_size ; (objective_function(x, y+grid_size) - objective_function(x, y))/ grid_size] ;

        %beta = (transpose(gamma) * linear_A * P) / (transpose(P) * A * P)
        beta = (transpose(gradient) * (gradient - gradient_old)) / (transpose((gradient - gradient_old)) * gradient_old) ;
        P = -gradient + beta * P ;

        iteration = iteration +1;
        x_all_iteration=[x_all_iteration x];
        y_all_iteration=[y_all_iteration y];
    end
    
    fprintf('--------------------------------------------\n')
    fprintf('the final (x, y) in nonlinear conjugate gradient method (HS): (%i, %i)\n', x, y);
    fprintf('the termination condition: %f\n', termination_condition)
    fprintf('the number of iteration: %i\n', iteration)
    toc
end

