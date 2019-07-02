function [ J ] = map_jacobian( a0, a1, mu, nu_vec)
%JACOBIAN_MAP Summary of this function goes here
%   Detailed explanation goes here

a = a0 + a1;

[Nmu, Nnu] = size(mu);
sh  = sinh(mu);
ch  = cosh(mu);
s = repmat( sin(nu_vec), Nmu, 1);
c = repmat( cos(nu_vec), Nmu, 1);


J = zeros(Nmu, Nnu, 3, 3);

J(:,:,1,1) = a.*ch.*s;
J(:,:,1,3) = a.*sh.*c;
J(:,:,2,1) = a.*sh.*c;
J(:,:,2,3) = -a.*ch.*s;
J(:,:,3,2) = a.*sh.*s;

end

