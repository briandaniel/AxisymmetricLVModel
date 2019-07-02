function [ Se ] = elastic_stress ( C_prol, E_prol, cvec, Q )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[ C_fib ] = tensor_rotate_prolate_to_fiber( Q, C_prol );
[ E_fib ] = tensor_rotate_prolate_to_fiber( Q, E_prol );


c1 = cvec(1);
c2 = cvec(2);
c3 = cvec(3);
c4 = cvec(4);

bff = 4*cvec(2) + 4*cvec(3);
bxx = 4*cvec(3);
bfx = 2*cvec(4) + 4*cvec(3);

eW = exp( bff.*E_fib(:,:,3,3).^2 + bxx.*(E_fib(:,:,2,2).^2 ...
          + E_fib(:,:,1,1).^2 + 2.*E_fib(:,:,1,2).^2) ...
          + bfx.*(2*E_fib(:,:,3,1).^2 + 2.*E_fib(:,:,3,2).^2 ) );

Se_fib = zeros(size(C_fib));     

Se_fib(:,:,1,1) = 2*c1*c3*eW.*(C_fib(:,:,1,1)-1);
Se_fib(:,:,2,2) = 2*c1*c3*eW.*(C_fib(:,:,2,2)-1);
Se_fib(:,:,3,3) = 2*c1*( c2 + c3 )*eW.*(C_fib(:,:,3,3)-1);

% Se_fib(:,:,1,2) = 4*c1*c3*eW.*C_fib(:,:,1,2);
% Se_fib(:,:,3,1) = 2*c1*( 2*c3 + c4 )*eW.*C_fib(:,:,3,1);
% Se_fib(:,:,3,2) = 2*c1*( 2*c3 + c4 )*eW.*C_fib(:,:,3,2);

Se_fib(:,:,1,2) = 2*c1*c3*eW.*C_fib(:,:,1,2);
Se_fib(:,:,3,1) = 1*c1*( 2*c3 + c4 )*eW.*C_fib(:,:,3,1);
Se_fib(:,:,3,2) = 1*c1*( 2*c3 + c4 )*eW.*C_fib(:,:,3,2);

Se_fib(:,:,2,1) = Se_fib(:,:,1,2);
Se_fib(:,:,1,3) = Se_fib(:,:,3,1);
Se_fib(:,:,2,3) = Se_fib(:,:,3,2);

% 
% eW = ( bff.*E_fib(:,:,3,3).^2 + bxx.*(E_fib(:,:,2,2).^2 ...
%           + E_fib(:,:,1,1).^2 + 2.*E_fib(:,:,1,2).^2) ...
%           + bfx.*(2*E_fib(:,:,3,1).^2 + 2.*E_fib(:,:,3,2).^2 ) );
% 
% Se_fib = zeros(size(C_fib));     
% 
% Se_fib(:,:,1,1) = c1*eW.*bxx.*E_fib(:,:,1,1);
% Se_fib(:,:,2,2) = c1*eW.*bxx.*E_fib(:,:,2,2);
% Se_fib(:,:,3,3) = c1*eW.*bff.*E_fib(:,:,3,3);
% Se_fib(:,:,1,2) = c1*eW.*bxx.*E_fib(:,:,1,2);
% Se_fib(:,:,3,1) = c1*eW.*bfx.*E_fib(:,:,3,1);
% Se_fib(:,:,3,2) = c1*eW.*bfx.*E_fib(:,:,3,2);
% 
% Se_fib(:,:,2,1) = Se_fib(:,:,1,2);
% Se_fib(:,:,1,3) = Se_fib(:,:,3,1);
% Se_fib(:,:,2,3) = Se_fib(:,:,3,2);

% Rotate from fiber coordinates to prolate spheroidal coordinates

Se = tensor_rotate_fiber_to_prolate(Q,Se_fib);

end








