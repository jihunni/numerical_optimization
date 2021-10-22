function result = objective_function_3(x,y)
    result = (1.5-x+x.*y).^2 + (2.25-x+x.*y.^2).^2 + (2.625-x+x.*y.^3).^2;
end