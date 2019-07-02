function [ Cinv ] = inverse_cauchy_green_deformation( C )
%INVERSE_CAUCHY_GREEN_DEFORMATION Summary of this function goes here
%   Detailed explanation goes here

Cinv = zeros(size(C));

detC = C(:,:,1,1).*C(:,:,2,2).*C(:,:,3,3) ...
       - C(:,:,1,1).*C(:,:,2,3).*C(:,:,3,2) ...
       - C(:,:,1,2).*C(:,:,2,1).*C(:,:,3,3);

Cinv(:,:,1,1) = 1./detC .*(  C(:,:,2,2).*C(:,:,3,3) - C(:,:,2,3).*C(:,:,3,2) );                            
Cinv(:,:,1,2) = 1./detC .*( - C(:,:,1,2).*C(:,:,3,3) );
Cinv(:,:,1,3) = 1./detC .*( C(:,:,1,2).*C(:,:,2,3) );

Cinv(:,:,2,1) = 1./detC .*( - C(:,:,3,3).*C(:,:,2,1) );
Cinv(:,:,2,2) = 1./detC .*( C(:,:,1,1).*C(:,:,3,3) );
Cinv(:,:,2,3) = 1./detC .*( - C(:,:,1,1).*C(:,:,2,3) );

Cinv(:,:,3,1) = 1./detC .*( C(:,:,2,1).*C(:,:,3,2) );
Cinv(:,:,3,2) = 1./detC .*( - C(:,:,1,1).*C(:,:,3,2) );
Cinv(:,:,3,3) = 1./detC .*(  C(:,:,1,1).*C(:,:,2,2) - C(:,:,1,2).*C(:,:,2,1) );

end


