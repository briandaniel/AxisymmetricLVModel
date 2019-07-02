function [ E, dE_da1, dE_da2, dE_da3 ] ...
             = green_strain_derivatives( F, dF_da1, dF_da2, dF_da3, C )
% green_strain_derivatives
% Calculates the green strain tensor { E = 1/2(C-I) } and derivatives in ai

dE_da1 = zeros(size(F));
dE_da2 = zeros(size(F));
dE_da3 = zeros(size(F));
E = zeros(size(F));

E(:,:,1,1) = .5*(C(:,:,1,1) - 1);
E(:,:,1,2) = .5*(C(:,:,1,2));
E(:,:,2,2) = .5*(C(:,:,2,2) - 1);
E(:,:,3,2) = .5*(C(:,:,3,2));
E(:,:,3,3) = .5*(C(:,:,3,3) - 1);

dE_da1(:,:,1,1) = F(:,:,1,1).*dF_da1(:,:,1,1);
dE_da1(:,:,3,3) = F(:,:,3,3).*dF_da1(:,:,3,3);
dE_da1(:,:,2,2) = F(:,:,1,2).*dF_da1(:,:,1,2) + F(:,:,2,2).*dF_da1(:,:,2,2) + F(:,:,3,2).*dF_da1(:,:,3,2);
dE_da1(:,:,1,2) = .5*( F(:,:,1,2).*dF_da1(:,:,1,1) + F(:,:,1,1).*dF_da1(:,:,1,2) );
dE_da1(:,:,3,2) = .5*( F(:,:,3,2).*dF_da1(:,:,3,3) + F(:,:,3,3).*dF_da1(:,:,3,2) );

dE_da2(:,:,1,1) = F(:,:,1,1).*dF_da2(:,:,1,1);
dE_da2(:,:,3,3) = F(:,:,3,3).*dF_da2(:,:,3,3);
dE_da2(:,:,2,2) = F(:,:,1,2).*dF_da2(:,:,1,2) + F(:,:,2,2).*dF_da2(:,:,2,2) + F(:,:,3,2).*dF_da2(:,:,3,2);
dE_da2(:,:,1,2) = .5*( F(:,:,1,2).*dF_da2(:,:,1,1) + F(:,:,1,1).*dF_da2(:,:,1,2) );
dE_da2(:,:,3,2) = .5*( F(:,:,3,2).*dF_da2(:,:,3,3) + F(:,:,3,3).*dF_da2(:,:,3,2) );

dE_da3(:,:,1,1) = F(:,:,1,1).*dF_da3(:,:,1,1);
dE_da3(:,:,3,3) = F(:,:,3,3).*dF_da3(:,:,3,3);
dE_da3(:,:,2,2) = F(:,:,1,2).*dF_da3(:,:,1,2) + F(:,:,2,2).*dF_da3(:,:,2,2) + F(:,:,3,2).*dF_da3(:,:,3,2);
dE_da3(:,:,1,2) = .5*( F(:,:,1,2).*dF_da3(:,:,1,1) + F(:,:,1,1).*dF_da3(:,:,1,2) );
dE_da3(:,:,3,2) = .5*( F(:,:,3,2).*dF_da3(:,:,3,3) + F(:,:,3,3).*dF_da3(:,:,3,2) );


% Fill in symmetric off-diagonals
E(:,:,2,1) = E(:,:,1,2);
E(:,:,2,3) = E(:,:,3,2);
dE_da1(:,:,2,1) = dE_da1(:,:,1,2);
dE_da1(:,:,2,3) = dE_da1(:,:,3,2);
dE_da2(:,:,2,1) = dE_da2(:,:,1,2);
dE_da2(:,:,2,3) = dE_da2(:,:,3,2);
dE_da3(:,:,2,1) = dE_da3(:,:,1,2);
dE_da3(:,:,2,3) = dE_da3(:,:,3,2);


end

