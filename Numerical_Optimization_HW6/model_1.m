function [value] = model_1(a,b,c,d, X,Y,Z)
%MODEL_1 Summary of this function goes here
%   Input : 
%            a : (float) model variable 
%            b : (float) model variable 
%            c : (float) model variable
%            d : (float) model variable
%            x : (float column vector, m*1) data
%            y : (float column vector, m*1) data
%            z : (float column vector, m*1) data
%   output : 
%            value : (float column vector, m*1) model value
%   Detailed explanation goes here
    value = [X, Y, Z, ones(length(X),1)] * [a ; b ; c ; d];
end

