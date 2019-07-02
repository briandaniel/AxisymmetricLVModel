function [eta1, eta2] = traction_integral( subel, param, nu_vec, a0, a1, a2  )
%TRACTION_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here


% Call backs
dmu_dnu  = subel.dmu_dnu ;
dmu_da1  = subel.dmu_da1 ;
dmu_da2 = subel.dmu_da2;
mu = subel.mu;
nu_vec = param.nu_vec; 
a0 = param.a0;


% Assign the inner wall value of \mu to use in integral calculation
mu_in = mu(1,:);

a = a0+a1;

dr_da1 = sin(nu_vec).*(sinh( mu_in ) + a*cosh(mu_in).*dmu_da1(1,:) );
dz_da1 = cos(nu_vec).*(cosh( mu_in ) + a*sinh(mu_in).*dmu_da1(1,:) );

dr_da2 = a*cosh(mu_in).*sin(nu_vec).*dmu_da2(1,:);
dz_da2 = a*sinh(mu_in).*cos(nu_vec).*dmu_da2(1,:);

nr = sqrt(2)*(cosh(mu_in).*sin(nu_vec) - cos(nu_vec).*sinh(mu_in).*dmu_dnu(1,:) ) ...
     ./sqrt( ( cosh(2*mu_in) - cos(2*nu_vec) ).*( 1 + dmu_dnu(1,:).^2 ));


nz = sqrt(2)*(sinh(mu_in).*cos(nu_vec) + cosh(mu_in).*sin(nu_vec).*dmu_dnu(1,:) ) ...
     ./sqrt(( cosh(2*mu_in) - cos(2*nu_vec) ).*( 1 + dmu_dnu(1,:).^2 ));

variable_changer = 1/sqrt(2)*sinh(mu_in).*sin(nu_vec).*sqrt( ( cosh(2*mu_in) - cos(2*nu_vec) ).*...
                        ( 1+ dmu_dnu(1,:).^2 ) );

integrand1 = (nr.*dr_da1 + nz.*dz_da1).*variable_changer;
integrand2 = (nr.*dr_da2 + nz.*dz_da2).*variable_changer;
 
eta1 = 2*pi*a^2*numeric_integral(nu_vec,integrand1);
eta2 = 2*pi*a^2*numeric_integral(nu_vec,integrand2);

end

