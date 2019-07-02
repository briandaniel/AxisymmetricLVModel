function [ wallVol ] = wall_volume_calc(  mu_in, mu_out, nu0, a)
%VOLUME_CALC Summary of this function goes here
%   Detailed explanation goes here


[x_in, y_in, z_in,x_out, y_out, z_out] = mu2xyz( mu_in, mu_out, nu0,0, a);

volout = pi*trapz(-z_out,x_out.^2);
volin = pi*trapz(-z_in,x_in.^2);

wallVol = volout-volin;

end

