function [x, y, x_all_iteration, y_all_iteration] = linear_conjugate_gradient_method_practical(objective_function, linear_A, linear_b, x, y, termination_condition)
%practical_conjugate_gradient_method Summary of this function goes here
% Multivariate optimization algorithm ; Netwon's method is a gradient-based
% optimization method. This program is implemented for 2 dimensional
% objective function f(x,y)
%   Detailed explanation goes here
%   Input : 
%            objective_function : (function) to minimize
%            linear_A : (float) 2x2 matrix, linear objective function in matrix
%            form
%            linear_b : (float) 2X1 matrix, linear objective function in matrix
%            form
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
    gamma = linear_A * [x ; y] - linear_b ; %2x1 matrix
    P = -gamma;

    %while(iteration < 5 | ...
    %        sqrt((x-x_all_iteration(end-1))^2 + (y-y_all_iteration(end-1)^2)) > termination_condition ...
    %      )
    while(transpose(gamma) * gamma >= termination_condition)
        %alpha = -(transpose(gamma) * P) / (transpose(P) * linear_A * P) ;
        alpha = (transpose(gamma) * gamma) / (transpose(P) * linear_A * P) ;
        x = x + alpha * P(1) ;
        y = y + alpha * P(2) ;
        gamma_old = gamma;
        gamma = linear_A * [x ; y] - linear_b;
        %beta = (transpose(gamma) * linear_A * P) / (transpose(P) * A * P)
        beta = (transpose(gamma) * gamma) / (transpose(gamma_old) * gamma_old) ;
        P = -gamma + beta * P ;

        iteration = iteration +1;
        x_all_iteration=[x_all_iteration x];
        y_all_iteration=[y_all_iteration y];
    end
    
    fprintf('--------------------------------------------\n')
    fprintf('the final (x, y) in linear conjugate gradient method (practical): (%i, %i)\n', x, y)
    fprintf('the termination condition: %f\n', termination_condition)
    fprintf('the number of iteration: %i\n', iteration)
    toc
end

