function [x, y, x_all_iteration, y_all_iteration] = quasi_netwons_method_SR1(objective_function, x, y, termination_condition)
%QUASI_NETWONS_METHOD Summary of this function goes here
% Multivariate optimization algorithm ; steepest_descent is a gradient-based
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
    grid_interval = 0.01;
    x_all_iteration=[x];
    y_all_iteration=[y];
    P = [5000 5000];
    objective_function_gradient = [0 ;0];
    B_inverse = eye(2) ; 
    
    while(iteration < 5 | ...
            sqrt((x-x_all_iteration(end-1))^2 + (y-y_all_iteration(end-1)^2)) > termination_condition ...
          )
        %Compute search direction P
        %% to update gradient on (x, y)
        objective_function_gradient_previous = objective_function_gradient;
        gradient_x = (objective_function(x+grid_interval, y) - objective_function(x,y))/grid_interval;
        gradient_y = (objective_function(x, y+grid_interval) - objective_function(x,y))/grid_interval;
        objective_function_gradient = [gradient_x; gradient_y];
        
        %% compute search direction P and search length alpha
        P = (-1) * B_inverse * objective_function_gradient ; % column vector
        alpha = exact_line_search(objective_function, x, y, P, objective_function_gradient);
        
        %% update the point
        x = x + alpha * P(1);
        y = y + alpha * P(2);        
        x_all_iteration = [x_all_iteration x];
        y_all_iteration = [y_all_iteration y];        

        % Compute matrix B^(-1) (approximation of hessian)
        I = eye(2); % identity matrix
        S = ([x_all_iteration(end) - x_all_iteration(end-1); y_all_iteration(end) - y_all_iteration(end-1)]) ; % column vector
        Y = objective_function_gradient - objective_function_gradient_previous ;
        
        B_inverse = B_inverse + (S-B_inverse*Y)*transpose(S-B_inverse*Y)/(transpose(S-B_inverse*Y)*Y);

        % output
        iteration = iteration + 1;
        fprintf('the current iteration in steepest desent: %i\n', iteration)
        fprintf('the current alpha in steepest desent: %i\n', alpha)
        fprintf('the current gradient in steepest desent: (%i, %i)\n', P(1), P(2))
        fprintf('the current x in steepest desent: %i\n', x)
        fprintf('the current y in steepest desent: %i\n', y)
        fprintf('\n')
    end
    
    fprintf('--------------------------------------------\n', iteration)
    fprintf('the final (x, y) in steepest desent: (%i, %i)\n', x, y)
    fprintf('the termination condition: %f\n', termination_condition)
    fprintf('the number of iteration: %i\n', iteration)
    toc
end

