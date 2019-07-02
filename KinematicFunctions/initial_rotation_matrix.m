function [ Q, omega ] ... 
    = initial_rotation_matrix( subel,param )
%ROTATION_MATRIX Summary of this function goes here
%   Detailed explanation goes here

% Call backs
sh0 = subel.sh0;
ch0 = subel.ch0;
s = subel.s;
c = subel.c;
Nmu = param.Nmu;
Nnu = param.Nnu;
psi_in_b0 = param.psi_in_b0;
psi_out_b0 = param.psi_out_b0;
muin0 = param.muin0;
muout0 = param.muout0;
muvec0 = subel.muvec0;

% Calculate the initial psi angle values as a function of initial mu values
psi_b0 = 1./(muout0 - muin0).*( psi_in_b0.*(muout0 - muvec0) ...
                                  + psi_out_b0.*(muvec0 - muin0) );
     
% Calculate angle parameter, omega
omega = 1./(tan(psi_b0).*tanh(muvec0));
omg = repmat( omega', 1, Nnu );

% Calculate psi
cpsi = sh0.*s.^2.*omg./sqrt( sh0.^2.*c.^2 + ch0.^2.*s.^2 ...
                                        + sh0.^2.*s.^4.*omg.^2 );                                    
psi = acos( cpsi );

Q = zeros(Nmu,Nnu,3,3);
Q(:,:,1,1) = 1;
Q(:,:,2,2) = cos(psi);
Q(:,:,2,3) = sin(psi);
Q(:,:,3,2) = -sin(psi);
Q(:,:,3,3) = cos(psi);




end


