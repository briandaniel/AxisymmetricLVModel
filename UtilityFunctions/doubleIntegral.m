function [ value ] = doubleIntegral( x, y, integrand )
%DOUBLEINTEGRAL Summary of this function goes here
%   x is the vector of x-values
%   y is the vector of y-values
%   The shape of the integrand should be have x increase across the
%   columns and y increases down the rows

integrand2 = zeros(size(x));

for k = 1:length(x);
    integrand2(k) = numeric_integral(y,integrand(k,:));
end

value = numeric_integral(x,integrand2);

end

