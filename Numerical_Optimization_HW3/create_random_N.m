function [X,Y] = create_random_N(N, range)
%CREATE_RANDOM_N Summary of this function goes here
%   Detailed explanation goes here
X = randi([-range range],1,N);
Y = randi([-range range],1,N);
end

