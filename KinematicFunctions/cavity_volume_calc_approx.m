function [ vol ] = cavity_volume_calc_approx(  mu_in, mu_out, nu0, a)
%VOLUME_CALC Summary of this function goes here
%   Detailed explanation goes here


[x_in, y_in, z_in,x_out, y_out, z_out] = mu2xyz( mu_in, mu_out, nu0,0, a);

vol = pi*trapz(-z_in,x_in.^2);

end

