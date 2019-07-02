function [ At] = activation_func( t, Ta, Tc, eps_fed, kd )
%ACTIVATION_FUNC Summary of this function goes here
%   Detailed explanation goes here

d = 1./(1+kd*eps_fed);

At = 0;

if t < Ta 
    At = ( sin(pi*t/Ta) ).^d;
end

end

