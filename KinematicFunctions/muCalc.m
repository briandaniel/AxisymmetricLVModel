function [ mu ] =  muCalc( a0, a1, a2, mu0, nu, mu_in0)
%MUCALC calculates the mu_in and mu_out given
%   a2, a, 


% Original formulation
    a = a0+a1;
    C0 = 1/3 - (cos(nu)).^2;   % Constant factor

    F = a2*a0^2*( cosh(mu_in0)*( 1/3*(cosh(mu_in0)).^2 - (cos(nu)).^2 ) - C0 );

    % Calculate mu_in    
    RHS = a0.^3 * ( cosh(mu0)*( 1/3*(cosh(mu0)).^2 - (cos(nu)).^2 ) - C0 ) -  F;   
    S = -(RHS/a^3 + C0);
    
    r3 = 2^(1/3);
    b = -(cos(nu)).^2;
    
    root3 = ( sqrt( 4*b.^3 + 9*S.^2 ) - 3*S ).^(1/3);
    mu = real( acosh( root3 ./ r3 - r3 .* b ./ root3 ) ); 

    % Formulation built to match mikes code
%     a = a0+a1;
%     c1 = 1/3*a.^3;
%     c2 = -a.^3.*cos(nu).^2;
%     
%     D0 = (cosh(mu_in0) + cosh(mu_in0).^2 + 1);
%     F = a2*(1-3/D0*cos(nu).^2);    
% %         F = a2*a0^2*( cosh(mu_in0)*( 1/3*(cosh(mu_in0)).^2 - (cos(nu)).^2 ) - C0 );
%     c3 = -a0.^3 * ( 1/3* cosh(mu0).^3 - cosh(mu0).*cos(nu).^2 + cos(nu).^2 - 1/3) + a.^3*(cos(nu).^2 - 1/3) + F;
%     box = ( sqrt(3)*sqrt(27*c1.^4.*c3.^2 + 4*c1.^3.*c2.^3) - 9*c1.^2.*c3).^(1/3);
%     root1 = box ./ (2.^(1/3) *3^(2/3)*c1) - (2/3).^(1/3).*c2./box;
% %     root2 = (1+1i*sqrt(3))*c2 ./(2^(2/3)*3^(1/3)*box) - (1-1i*sqrt(3))*box./(2^(4/3).*3^(2/3)*c1);
% %     root3 = (1-1i*sqrt(3))*c2 ./(2^(2/3)*3^(1/3)*box) - (1+1i*sqrt(3))*box./(2^(4/3).*3^(2/3)*c1);
%     pos1 = acosh(root1);
% %     pos2 = acosh(root2);
% %     pos3 = acosh(root3);
%     
%     mu = real(pos1);

end
    

