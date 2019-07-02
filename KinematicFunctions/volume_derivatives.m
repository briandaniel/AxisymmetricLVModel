function [ dV_da1, dV_da2 ] = volume_derivatives( subel,param)
%VOLUME_DERIVATIVES Summary of this function goes here
%   Detailed explanation goes here


% Call backs
a = subel.a;
dmu_dnu  = subel.dmu_dnu ;
dmu_da1  = subel.dmu_da1 ;
dmu_da2 = subel.dmu_da2;
ddmu_dnuda1 = subel.ddmu_dnuda1;
ddmu_dnuda2  = subel.ddmu_dnuda2 ;
mu = subel.mu;
nu_vec = param.nu_vec; 


% values only taken on inner wall
mu = mu(1,:);
nu = nu_vec;
dmu_dnu = dmu_dnu(1,:);
ddmu_dnuda1 = ddmu_dnuda1(1,:);
ddmu_dnuda2 = ddmu_dnuda2(1,:);
dmu_da2 = dmu_da2(1,:);
dmu_da1 = dmu_da1(1,:);

s = sin(nu);
s_sq = s.^2;
sh = sinh(mu);
c = cos(nu);
sh_sq = sh.^2;

integrand1 = a^2/2*s_sq.*sh.*(-3.*s.*sinh(2*mu)+2*c.*sh_sq.*(3*dmu_dnu + a.*ddmu_dnuda1) ...
                - a*(s+3*cosh(2*mu).*s - 3.*c.*dmu_dnu.*sinh(2*mu) ).*dmu_da1);

integrand2 = a^3/2*s_sq.*sh.*(-2.*c.*sh_sq.*ddmu_dnuda2 + ...
            (s + 3*cosh(2*mu).*s - 3*dmu_dnu.*c.*sinh(2*mu)).*dmu_da2);
        
dV_da1 = pi* numeric_integral(nu_vec,integrand1);
dV_da2 = pi* numeric_integral(nu_vec,integrand2);


end

