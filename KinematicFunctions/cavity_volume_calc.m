function [ vol ] = cavity_volume_calc(  mu_in, nu_vec, dmu_dnu_in  , a)
%VOLUME_CALC Summary of this function goes here
%   Detailed explanation goes here

integrand = sinh(mu_in).^2.*sin(nu_vec).^2.*(sinh(mu_in).*cos(nu_vec).*dmu_dnu_in ...
                - cosh(mu_in).*sin(nu_vec));

vol = pi*a.^3.*numeric_integral(-nu_vec,integrand);

end

