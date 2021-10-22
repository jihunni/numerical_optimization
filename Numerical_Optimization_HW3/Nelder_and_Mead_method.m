function [x,y, x_all_iteration, y_all_iteration] = Nelder_and_Mead_method(x, y, objective_function, alpha, beta, gamma, termination_condition)
%NELDER_AND_MEAD_METHOD Summary of this function goes here
% INPUTS
% x
% y
% objective_function
% alpha
% beta
% gamma
% termination_condition : area of simplex

%   Detailed explanation goes here
% To measure the performance
tic
iteration = 0;

x_all_iteration=[];
y_all_iteration=[];
%objective_function = @(x,y)objective_function(x, y);
%check the parameter
if (alpha > 0 && beta > 1 && 0 < gamma && gamma < 1)
else
    fprintf("check parameter again\n")
end

while(polyarea(x, y) > termination_condition)
    % sort and set up for Reflection
    

    function_value = objective_function(x,y);
    [function_value, index_sorted]=sort(function_value);
    x = x(index_sorted);
    y = y(index_sorted);
    
    % Reflection
    centroid_x = sum(x(1:end-1))/length(x(1:end-1));
    centroid_y = sum(y(1:end-1))/length(y(1:end-1));
    x_r = centroid_x + alpha * (centroid_x-x(end));
    y_r = centroid_y + alpha * (centroid_y-y(end));
    function_value_r = objective_function(x_r, y_r);
    
    if(function_value(1) <= function_value_r && function_value_r <= function_value(end-1))
        x(end) = x_r;
        y(end) = y_r;
        function_value(end) = function_value_r;
    
    elseif(function_value_r >= function_value(end-1))
        %contraction
        if (function_value_r < function_value(end))
            x_c = centroid_x + gamma * (x_r - centroid_x);
            y_c = centroid_y + gamma * (y_r - centroid_x);
        else
            x_c = centroid_x + gamma * (x(end) - centroid_x);
            y_c = centroid_y + gamma * (y(end) - centroid_x);
        end

        function_value_c = objective_function(x_c, y_c);

        if (function_value_c < min(function_value_r, function_value(end)))
            x(end) = x_c;
            y(end) = y_c;
        else
            x_1 = x(1);
            x = (x+x_1)./2;
        end

    elseif(function_value_r <= function_value(1))
        %expansion
        x_e = centroid_x + beta * (x_r - centroid_x);
        y_e = centroid_y + beta * (x_r - centroid_y);
        function_value_e = objective_function(x_e, y_e);
    
        if (function_value_e <= function_value_r)
            x(end) = x_e;
            y(end) = y_e;
        else
            x(end) = x_r;
            y(end) = y_r;
        end
    else
        fprintf('ERROR!\n');
    end
    
    iteration = iteration + 1;

    x_all_iteration = [x_all_iteration x];
    y_all_iteration = [y_all_iteration y];
end %while loop
fprintf('the termination condition: %f\n', termination_condition)
fprintf('the number of iteration: %i\n', iteration)
toc
end

