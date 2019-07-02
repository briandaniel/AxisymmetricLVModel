function [ Sf_ai_prol ] = active_fiber_stress( At, kav, kp, m, ...
                        ell, ell0, eps_fed, dell_dai, Q, F )
%ACTIVE_FIBER_STRESS Summary of this function goes here
%   Detailed explanation goes here

eps_f = 1/2*(ell.^2./ell0.^2 - 1 );

% active fiber stress component in the ff direction (All other directions
% are zero)
sigmaF = At.*kav.*(1+m.*eps_f).*(1+kp.*eps_fed).*ell./ell0.^2.*dell_dai;

% Convert to Piola-Kirchoff stress tensor
SF = ell0.^2./ell.^2 .*sigmaF;

% Rotate from fiber coordinates back to prolate spheroidal coordaintes
[Nmu,Nnu] = size(ell);
Sf_ai_fib = zeros(Nmu,Nnu,3,3);
Sf_ai_fib(:,:,3,3) = SF;

Sf_ai_prol = tensor_rotate_fiber_to_prolate(Q,Sf_ai_fib);


end

