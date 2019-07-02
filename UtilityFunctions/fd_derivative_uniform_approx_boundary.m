function [ dydx ] = fd_derivative_uniform_approx_boundary( y, dx, derivative_number )
% Calculates the n^th derivative (stored in derivative_number) of y using 
% a uniform grid with spacing dx. The stencil uses centered finite
% difference except at the boundaries which are approximated by the
% corresponding left and right hand rules.

% Currently computes only second order

dydx = zeros(size(y));

if derivative_number == 1
    dydx(2:end-1) = ( y(3:end) - y(1:end-2) ) / (2*dx);
    dydx(1) = 1/dx* ( -3/2*y(1) + 2*y(2) -1/2*y(3) );
    dydx(end) = 1/dx* ( 3/2*y(end) - 2*y(end-1) + 1/2*y(end-2) );
end


end

