function [ ell, ell0, dell_da1, dell_da2, dell_da3 ] ...
    = ell_derivative( subel, param, a0, a1, a2, a3, omega )
%ELL_DERIVATIVES Summary of this function goes here
%   Detailed explanation goes here

% Call backs
a = subel.a;
nu = subel.nu;
sh_sq = subel.sh_sq;
ch_sq = subel.ch_sq;
sh0 = subel.sh0;
ch0 = subel.ch0;
s = subel.s;
c_sq = subel.c_sq;
s_sq = subel.s_sq;
dmu_dnu  = subel.dmu_dnu ;
dmu_da1  = subel.dmu_da1 ;
dmu_da2 = subel.dmu_da2;
ddmu_dnuda1 = subel.ddmu_dnuda1;
ddmu_dnuda2  = subel.ddmu_dnuda2 ;
mu = subel.mu;
Nnu = param.Nnu;
a0 = param.a0;


omg = repmat(omega',1,Nnu);

ell = a*sqrt( (1+dmu_dnu.^2).*(ch_sq.*s_sq + sh_sq.*c_sq) ...
                                        + (omg+a3).^2.*sh_sq.*s.^4 );

ell0 = a0 * sqrt( ch0.^2.*s_sq + sh0.^2.*c_sq + omg.^2.*sh0.^2.*s.^4);

denom = sqrt(-2*(cos(2*nu) - cosh(2*mu) ).*(1 + dmu_dnu.^2 ) ...
                + 4*(a3 + omg).^2.*s.^4.*sh_sq );

dell_da1 = ( - cos(2*nu) + cosh(2*mu) + 2.*(a3+omg).^2.*s.^4.*sh_sq - ...
    (cos(2*nu) - cosh(2*mu)).*dmu_dnu.*(dmu_dnu+a.*ddmu_dnuda1) + ...
    a.*(1 + dmu_dnu.^2 + (a3+omg).^2.*s.^4).*sinh(2*mu).*dmu_da1 )./denom;

dell_da2 = a.*( (-cos(2*nu) + cosh(2*mu) ).*dmu_dnu.*ddmu_dnuda2 ...
    + 1/8.*(8 + 3*(a3+omg).^2 + (a3 + omg).^2.*(-4.*cos(2*nu)+cos(4*nu)) ...
    + 8.*dmu_dnu.^2 ) .*sinh(2*mu).*dmu_da2  )./denom;

dell_da3 = 2*a*(a3+omg).*s.^4.*sh_sq./denom;

end

