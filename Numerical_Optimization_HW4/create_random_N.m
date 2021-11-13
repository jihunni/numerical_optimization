function [X,Y] = create_random_N(N, range)
%CREATE_RANDOM_N Summary of this function goes here
%   Detailed explanation goes here
%   INPUT : 
%           N : integer, the number of pair (x, y)
%           range : float, the random range
%   OUTPUT : 
%           X :  N-dimensional float vector
%           Y :  N-dimensional float vector
X = randi([-range range],1,N);
Y = randi([-range range],1,N);
end

