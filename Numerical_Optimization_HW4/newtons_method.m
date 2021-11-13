function [x, y, x_all_iteration, y_all_iteration] = newtons_method(objective_function, x, y, termination_condition)
%Newton's method Summary of this function goes here
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
    grid_interval = 0.01;
    x_all_iteration=[x];
    y_all_iteration=[y];
    P = [5000 5000];
    objective_function_gradient = [5000 5000];
    
    while(iteration < 5 | ...
            sqrt((x-x_all_iteration(end-1))^2 + (y-y_all_iteration(end-1)^2)) > termination_condition ...
          )
        % objective_function_gradient * transpose(objective_function_gradient) > termination_condition


        %% grid
        %X_grid = -abs(abs(x)+grid_interval*2):grid_interval:abs(abs(x)+grid_interval*2) ;
        %Y_grid = -abs(abs(y)+grid_interval*2):grid_interval:abs(abs(y)+grid_interval*2) ;
        %Y_grid = transpose(Y_grid);
        %objective_function_grid = objective_function(X_grid, Y_grid);
        %[fx, fy] = gradient(objective_function_grid, grid_interval);

        %% calculate the gradient on (x, y)
        %t = (X_grid == x) & (Y_grid == y);
        %indt = find(t);
        %objective_function_gradient = [fx(indt) fy(indt)] ;
        
        %if (indt == 0)
        %    fprintf("no index finding on gradient grid\n")
        %    %break
        %end

        %% to approximate hessian on (x, y) - next point
        %X_grid = -abs(abs(x)+grid_interval*2):grid_interval:abs(abs(x)+grid_interval*2) ;
        %Y_grid = -abs(abs(y)+grid_interval*2):grid_interval:abs(abs(y)+grid_interval*2) ;
        %Y_grid = transpose(Y_grid);

        %t_x = (X_grid == x+grid_interval) & (Y_grid == y);
        %indt_x = find(t_x);
        %t_y = (X_grid == x) & (Y_grid == y+grid_interval);
        %indt_y = find(t_y);

        %if (indt_x == 0 | indt_y == 0)
        %    fprintf("no index finding on hessian grid\n")
        %    break
        %end
        %hessian = [(fx(indt_x) - fx(indt)) / grid_interval ...
        %            (fx(indt_y) - fx(indt)) / grid_interval; ...
        %            (fx(indt_x) - fx(indt)) / grid_interval ...
        %            (fy(indt_y) - fy(indt)) / grid_interval] ;
        %hessian_inverse = inv(hessian); % 2*2 matrix
        %clear t; clear indt; clear t_x; clear indt_x; clear t_y; clear indt_y;

        % to approximate gradient
        %% gradient on (x, y)
        gradient_x0 = (objective_function(x+grid_interval, y) - objective_function(x,y))/grid_interval;
        gradient_y0 = (objective_function(x, y+grid_interval) - objective_function(x,y))/grid_interval;
        objective_function_gradient0 = [gradient_x0, gradient_y0];
        %% gradient on (x+grid_interval, y)
        gradient_x_dx = (objective_function(x+grid_interval*2, y) - objective_function(x+grid_interval,y))/grid_interval;
        gradient_y_dx = (objective_function(x+grid_interval, y+grid_interval) - objective_function(x+grid_interval,y))/grid_interval;
        objective_function_gradient_dx = [gradient_x_dx, gradient_y_dx];
        %% gradient on (x, y+grid_interval)
        gradient_x_dy = (objective_function(x+grid_interval, y+grid_interval) - objective_function(x,y+grid_interval))/grid_interval;
        gradient_y_dy = (objective_function(x, y+grid_interval*2) - objective_function(x,y+grid_interval))/grid_interval;
        objective_function_gradient_dx = [gradient_x_dy, gradient_y_dy];
        % to approximate hessian
        hessian = [(gradient_x_dx-gradient_x0)/grid_interval ...
                    (gradient_x_dy-gradient_y0)/grid_interval; ...
                    (gradient_y_dx-gradient_x0)/grid_interval ...
                    (gradient_y_dy-gradient_y0)/grid_interval];
        hessian_inverse = inv(hessian) ;
        % iteration for gradient-based method
        %gradient_direction = - objective_function_gradient0 ./ (objective_function_gradient0 * transpose(objective_function_gradient0)) ;
        %gradient_direction = transpose(gradient_direction) ; % to make it a column vector
        P = (-1) * hessian_inverse * transpose(objective_function_gradient0) ; % column vector
        alpha = inexact_line_search(objective_function, x, y, transpose(P), objective_function_gradient0);
       
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