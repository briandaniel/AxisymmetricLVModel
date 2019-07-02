function [ x_in, y_in, z_in, x_out, y_out, z_out] = ...
                        mu2xyz( mu_in, mu_out, nu0, phi, a)
%MU2XYZ Summary of this function goes here
%   Detailed explanation goes here

    % Calculate initial coordinates of inner wall
    x_in = a*sinh(mu_in).*sin(nu0).*cos(phi);
    y_in = a*sinh(mu_in).*sin(nu0).*sin(phi);
    z_in = a*cosh(mu_in).*cos(nu0);

    % Calculate initial coordinates of outer wall
    x_out = a*sinh(mu_out).*sin(nu0).*cos(phi);
    y_out = a*sinh(mu_out).*sin(nu0).*sin(phi);
    z_out = a*cosh(mu_out).*cos(nu0);

end

