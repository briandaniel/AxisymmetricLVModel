function [ integral ] = numeric_integral(x,f)
% THIS INTEGRATOR ASSUMES CONSTANT STEP SIZE!
%  x vector(or matrix) is the domain  
%  f gives the function values on x
%  dim gives which dimension to integrate across.

% integral = trapz(x,f)

dx = x(2)-x(1);

integral = dx/48*( 17*f(1) + 59*f(2) + 43*f(3) + 49*f(4) + 48*sum(f(5:end-4))...
                  +49*f(end-3) + 43*f(end-2) + 59*f(end-1) + 17*f(end) );
              
              
              
end

