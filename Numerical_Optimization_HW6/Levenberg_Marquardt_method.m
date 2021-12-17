function [a, b, c, d, x, y, z, x_all_iteration, y_all_iteration] = Levenberg_Marquardt_method(model, X, Y, Z, Val, termination_condition)
%GAUSS_NEWTON_METHOD Summary of this function goes here
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
    a = random('Normal',-10,10);
    b = random('Normal',-10,10);
    c = random('Normal',-10,10);
    d = random('Normal',-10,10);
    a_all_iteration=[a];
    b_all_iteration=[b];
    c_all_iteration=[c];
    d_all_iteration=[d];
    
    P = [5000 5000];
    lambda=random('Normal',0,10);

    while((iteration < 2 || ...
            sqrt((a-a_all_iteration(end-1))^2 + (b-b_all_iteration(end-1))^2 + (c-c_all_iteration(end-1))^2 + (d-d_all_iteration(end-1))^2) > termination_condition) & ...
          iteration < 100)
    % maximum iteraiton is 100
        % Compute Jacobian, gradient, hessian
        R_function = @(a,b,c,d)model(a, b, c, d, X,Y,Z) - Val ; % residual vector
        R = R_function(a,b,c,d);
        %R = model_in(a, b, c, d, X,Y,Z) - Val ; % residual vector
        %J = jacobian(R, [a, b, c, d]);%Jacobian
        J = [(R_function(a+grid_interval,b,c,d) - R_function(a,b,c,d))./grid_interval, (R_function(a,b+grid_interval,c,d) - R_function(a,b,c,d))./grid_interval, (R_function(a,b,c+grid_interval,d) - R_function(a,b,c,d))./grid_interval, (R_function(a,b,c,d+grid_interval) - R_function(a,b,c,d))./grid_interval;];
        gradient = transpose(J) * R;
        hessian = transpose(J) * J;

        %Compute search direction P(a,b,c,d)
        cal = transpose(J)*J;
        size_cal = size(cal);

        P = (-1) * inv(cal+lambda*eye(size_cal(1))) * transpose(J) * R ; % column vector
        
        % temporarily update the pointtoc
        a = a + P(1);
        b = b + P(2);        
        c = c + P(3);        
        d = d + P(4);        
        
        % Adjust lambda accordingly
        if sum(R_function(a,b,c,d).^2) > sum(R_function(a_all_iteration(end),a_all_iteration(end),c_all_iteration(end),d_all_iteration(end)).^2) 
        %if model(a,b,c,d, X, Y, Z) > model(a_all_iteration(end),b_all_iteration(end),c_all_iteration(end),d_all_iteration(end), X, Y, Z)
            % if not satisfy descent condition
            while sum(R_function(a,b,c,d).^2) >= sum(R_function(a_all_iteration(end),a_all_iteration(end),c_all_iteration(end),d_all_iteration(end)).^2) && lambda < 100
            % until to find the descendent condition
            %while model(a,b,c,d, X, Y, Z) >= model(a_all_iteration(end),b_all_iteration(end),c_all_iteration(end),d_all_iteration(end), X, Y, Z) 
                lambda = lambda * 10; % lambda
                P = (-1) * inv(cal+lambda*eye(size_cal(1))) * transpose(J) * R ; % column vector
                a = a_all_iteration(end) + P(1);
                b = b_all_iteration(end) + P(2);        
                c = c_all_iteration(end) + P(3);        
                d = d_all_iteration(end) + P(4);  
            end
        else 
            % if satisfy descent condition
            while sum(R_function(a,b,c,d).^2) < sum(R_function(a_all_iteration(end),a_all_iteration(end),c_all_iteration(end),d_all_iteration(end)).^2) && lambda > 0.1
                %until to the last descedent condition
                %while model(a,b,c,d, X, Y, Z) < model(a_all_iteration(end),b_all_iteration(end),c_all_iteration(end),d_all_iteration(end), X, Y, Z) 
                lambda = lambda * 0.1; % lambda
                P = (-1) * inv(cal+lambda*eye(size_cal(1))) * transpose(J) * R ; % column vector
                a = a_all_iteration(end) + P(1);
                b = b_all_iteration(end) + P(2);        
                c = c_all_iteration(end) + P(3);        
                d = d_all_iteration(end) + P(4);
            end
            % go back to previous step (last desending step)
            lambda = lambda * 10 ;
            P = (-1) * inv(cal+lambda*eye(size_cal(1))) * transpose(J) * R ; % column vector
            a = a_all_iteration(end) + P(1);
            b = b_all_iteration(end) + P(2);        
            c = c_all_iteration(end) + P(3);        
            d = d_all_iteration(end) + P(4);
        end
    
        a_all_iteration = [a_all_iteration a];
        b_all_iteration = [b_all_iteration b];
        c_all_iteration = [c_all_iteration c];
        d_all_iteration = [d_all_iteration d];
        
        % output
        iteration = iteration + 1;
        fprintf('the current iteration in steepest desent: %i\n', iteration)
        fprintf('the current gradient in steepest desent: (%i, %i, %i, %i)\n', P(1), P(2), P(3), P(4))
        fprintf('the current (a, b, c, d) in steepest desent: (%i, %i, %i, %i)\n', a, b, c, d)
        fprintf('the current lambda in steepest desent: (%i)\n', lambda)
        fprintf('\n')
    end
    
    cost = sum(R_function(a,b,c,d).^2) ;

    fprintf('--------------------------------------------\n', iteration)
    fprintf('the final (a, b, c, d) in steepest desent: (%i, %i, %i, %i)\n', a, b, c, d)
    fprintf('the final lambda in steepest desent: (%i)\n', lambda)
    fprintf('the final cost in steepest desent: %i\n', cost)
    fprintf('the termination condition: %f\n', termination_condition)
    fprintf('the number of iteration: %i\n', iteration)
    toc
end

