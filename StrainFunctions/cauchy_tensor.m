function [ C ] ...
                = cauchy_tensor( F )
% cauchy_tensor.m
% Calculates the cauchy deformation gradient tensor C = F^T*F
% Note that C21, C32 are nonzero, but the tensor is symmetric so these are
% equal to C12, C23, respectively.

C = zeros(size(F));

C(:,:,1,1) = F(:,:,1,1).^2;
C(:,:,1,2) = F(:,:,1,1).*F(:,:,1,2);
C(:,:,3,2)= F(:,:,3,3).*F(:,:,3,2);
C(:,:,3,3)= F(:,:,3,3).^2;
C(:,:,2,2)= F(:,:,1,2).^2 + F(:,:,2,2).^2 + F(:,:,3,2).^2;

C(:,:,2,3) = C(:,:,3,2);
C(:,:,2,1) = C(:,:,1,2);

end


