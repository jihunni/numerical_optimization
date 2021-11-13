function [x, y, x_all_iteration, y_all_iteration] = steepest_descent(objective_function, x, y, termination_condition)
%steepest_descent Summary of this function goes here
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
    grid_interval = 0.1;
    x_all_iteration=[x];iteration = 0;
    % initialization
    grid_interval = 0.1;
    x_all_iteration=[x];
    y_all_iteration=[y];
    P = [3 ; 3];
    
    
    while(iteration < 5 | sqrt((x-x_all_iteration(end-1))^2 + (y-y_all_iteration(end-1)^2)) > termination_condition) 
        %    sqrt(P(1)^2 + P(2)^2) > 0.01)
        % to find gradient
        %% grid
        %X_grid = -abs(x+grid_interval):grid_interval:abs(x+grid_interval) ;
        %Y_grid = -abs(y+grid_interval):grid_interval:abs(y+grid_interval) ;
        %Y_grid = transpose(Y_grid);
        %objective_function_grid = objective_function(X_grid, Y_grid);
        %[fx, fy] = gradient(objective_function_grid, grid_interval);
        %% calculate the gradient on (x, y)
        %t = (x == x) & (y == y);
        %indt = find(t);
        %objective_function_gradient = [fx(indt) fy(indt)] ;
        
        %% to approximate the gradient on (x, y)
        gradient_x0 = (double(objective_function(x+grid_interval, y)) - double(objective_function(x,y)))/grid_interval;
        gradient_y0 = (objective_function(x, y+grid_interval) - objective_function(x,y))/grid_interval;
        objective_function_gradient = [gradient_x0; gradient_y0];
      
        if (transpose(objective_function_gradient)*objective_function_gradient < 0.001)
            fprintf('break\n')
            break
        end

        % iteration for gradient-based method
        P = - objective_function_gradient ;
        alpha = exact_line_search(objective_function, x, y, P, objective_function_gradient);
        x = x + alpha * P(1);
        y = y + alpha * P(2);
            
        % output
        iteration = iteration + 1;
        x_all_iteration = [x_all_iteration x];
        y_all_iteration = [y_all_iteration y];

        
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