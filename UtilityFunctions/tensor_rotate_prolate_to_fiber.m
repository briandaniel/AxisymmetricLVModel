function [ Trot ] = tensor_rotate_prolate_to_fiber( Q, T )
% This file is model specific and assumes Q has the form
% Q = [ 1     0         0    ]
%     [ 0 cos(psi)  sin(psi) ]
%     [ 0 -sin(psi) cos(psi) ]
% 
% This file rotates 
%   FROM prolate coordinates 
%   TO fiber coordinates

cpsi = Q(:,:,2,2);
spsi = Q(:,:,2,3);

s2psi = 2.*spsi.*cpsi;
c2psi = cpsi.^2 - spsi.^2;

%%%%%%%%%%%%%%%%%%%%%%%%% signs adjusted to match mikes code throughout %%%%%%%%%%%%%%%%%%%%%%%%
% (this is just the switching of the sign between sin functions in the
% rotation matrix)

Trot(:,:,1,1) = T(:,:,1,1);
Trot(:,:,1,2) = T(:,:,1,2).*cpsi - T(:,:,1,3).*spsi;
Trot(:,:,1,3) = T(:,:,1,3).*cpsi + T(:,:,1,2).*spsi;

Trot(:,:,2,1) = Trot(:,:,1,2) ;

Trot(:,:,2,2) = T(:,:,2,2).*cpsi.^2 + (T(:,:,3,3).*spsi - 2*T(:,:,2,3).*cpsi).*spsi;
                
Trot(:,:,2,3) = T(:,:,2,3).*c2psi + (T(:,:,2,2) -  T(:,:,3,3)).*cpsi.*spsi;
                
Trot(:,:,3,1) = Trot(:,:,1,3);
Trot(:,:,3,2) = Trot(:,:,2,3);

Trot(:,:,3,3) = T(:,:,3,3).*cpsi.^2 + T(:,:,2,2).*spsi.^2 + s2psi.*T(:,:,2,3);
                    
 
% % VESTIGIAL CODE: MATRIX ROTATION TEST                
% [Nmu,Nnu,n,m] = size(Q);
% Qmtrx = normalMatrix(Q); 
% Tmtrx = normalMatrix(T);
% Trotmtrx = normalMatrix(Trot);
% 
% Temp = zeros(size(Qmtrx));
% 
% errors = zeros(Nmu,Nnu);
% for k = 1:Nmu, 
%     for j = 1:Nnu, 
%         Temp(:,:,k,j) = Qmtrx(:,:,k,j)'*Tmtrx(:,:,k,j)*Qmtrx(:,:,k,j);
%         
%         errors(k,j) = norm(Temp(:,:,k,j)-Trotmtrx(:,:,k,j));
%     end; 
% end; 
%          
% figure
% imagesc(errors)
% colorbar

end

