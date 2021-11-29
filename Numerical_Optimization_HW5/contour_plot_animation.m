function [] = contour_plot_animation(objective_function, x_iteration,y_iteration)
%CONTOUR_PLOT_ANIMATION Summary of this function goes here
%   Detailed explanation goes here
%   INPUT : 
%           objective_function : (function)
%           x_iteratione : (float list) the result of optimization
%           algorithm
%           y_iteratione : (float list) the result of optimization
%           algorithm
%   OUTPUT : 
%           Plot window
% grid
    x_grid = linspace(min(x_iteration) * 1.2,max(x_iteration) * 1.2);
    y_grid = linspace(min(y_iteration) * 1.2,max(y_iteration) * 1.2);
    [X, Y] = meshgrid(x_grid,y_grid);
    Z = objective_function(X,Y);
    
    % for animation plot
    for k = 1:(length(x_iteration)-1)
        % contour plot
        contour(X,Y,Z,10, 'ShowText','on')
        hold on
        
        % marker plot
        plot([x_iteration(1:k) x_iteration(1:k+1)], [y_iteration(1:k) y_iteration(1:k+1)], 'xr')
    
        % graph properties
        axis([min(x_iteration)*1.2 max(x_iteration)*1.2 min(y_iteration)*1.2 max(y_iteration)*1.2])
        grid on
        xlabel('x')
        ylabel('y')
        %legend('cos(t) marker', 'cos(t)')
        pause(1)
        
        if k ~= length(x_iteration)-1
            clf
        end
    end

end

