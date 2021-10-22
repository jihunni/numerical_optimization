function [x,y, x_all_iteration, y_all_iteration] = Powells_method(x_0, y_0, objective_function, termination_condition)
%POWELLS_METHOD Summary of this function goes here
%INPUT
%initial_pocint is a vector of input. (e.g.) (1,2,3)
%   Detailed explanation goes here
% (x, y) is p_0
% to measure the performacne
iteration=0;
tic

% initialization
direction_vector_1 = [1,0]; %vector u
direction_vector_2 = [0,1]; 
step_length=10000000000; % big number to make while loop works for first iteration
x = x_0;
y = y_0 ;
p_1 = [0,0];
p_2 = [0,0];

x_all_iteration=[];
y_all_iteration=[];
%objective_function(1,1)

while(step_length > termination_condition)
    p_0 = [x, y];
    
    % find direction vector
    serach_function = @(step_length)objective_function(p_0(1) + step_length*direction_vector_1(1), p_0(2) + step_length*direction_vector_1(2));
    step_length = fminsearch(serach_function, rand(1)*randi([0,10],1,1));
    %p_1 = p_0 + step_length .* direction_vector_1;
    p_1(1) = p_0(1) + step_length * direction_vector_1(1);
    p_1(2) = p_0(2) + step_length * direction_vector_1(2);
    


    serach_function = @(step_length)objective_function(p_1(1) + step_length*direction_vector_2(1), p_1(2) + step_length*direction_vector_2(2));
    step_length = fminsearch(serach_function, rand(1)*randi([0,10],1,1));
    p_2 = p_1 + step_length * direction_vector_2;
    
    direction_vector_1 = direction_vector_2;
    direction_vector_2(1) = p_2(1) - p_0(1);
    direction_vector_2(2) = p_2(2) - p_0(2);

    serach_function = @(step_length)objective_function(p_0(1) + direction_vector_2(1) * step_length, p_0(2) + direction_vector_2(2) * step_length);
    step_length = fminsearch(serach_function, rand(1)*randi([0,10],1,1));
    
    x = p_0(1) + direction_vector_2(1) * step_length;
    y = p_0(2) + direction_vector_2(2) * step_length; 

    
    iteration = iteration+1;

    %for animation graph
    x_all_iteration = [x_all_iteration x];
    y_all_iteration = [y_all_iteration y];

end %while loop
fprintf('the termination condition: %f\n', termination_condition)
fprintf('the number of iteration: %i\n', iteration)
toc
end

