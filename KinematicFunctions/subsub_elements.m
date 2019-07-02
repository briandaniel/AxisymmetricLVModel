function [ subel] ...
    = subsub_elements( mu, nu_vec, a0, a1, a2, mu0_vec, muin0 )
%SUBSUB_ELEMENTS Summary of this function goes here
%   Detailed explanation goes here

subel = struct();

a = a0 + a1;

[Nmu, Nnu] = size(mu);
nu =  repmat( nu_vec, Nmu, 1);

% Calculate basic elements to be used in further calculations
sh  = sinh(mu);
ch  = cosh(mu);
cth = coth(mu);
cch = csch(mu);
sh_sq = sh.^2;
ch_sq = ch.^2;
chin0 = cosh(muin0);
sh0  = repmat( sinh(mu0_vec)', 1, Nnu );
ch0  = repmat( cosh(mu0_vec)', 1, Nnu );

s = repmat( sin(nu_vec), Nmu, 1);
c = repmat( cos(nu_vec), Nmu, 1);
c_sq = c.^2;
s_sq = s.^2;
cch_sq = cch.^2;

sum_sq = sh_sq + s_sq;
sum0_sq = sh0.^2 + s_sq;
root_sum_sq = sqrt(sum_sq);
root_sum0_sq = sqrt(sum0_sq);

% Calculate subelements
% dmu_dmu0 has been checked
dmu_dmu0 = a0^3 * sh0 .* sum0_sq ...
                    ./ ( a^3 .* sh.* sum_sq ) ;     

% dmu_dnu has been checked
dmu_dnu = 2*c.*s.*(-a^3*(ch-1) + a0^3*(ch0 - 1) - a2*a0^2*(chin0 - 1) )...
                    ./ (a^3 * sh .* sum_sq );  
              
% dmu_da1 has been checked
dmu_da1 = (1 - ch.^3 + (3*ch - 3).*c_sq ) ...
                    ./(a*sh.*(ch_sq - c_sq) );
                
% dmu_da2 has been checked
dmu_da2 = a0^2*(1 - chin0^3 + (3*chin0 - 3)*c_sq )...
                    ./(3*a^3*sh.*(ch_sq - c_sq) );
                
% ddmu_dnuda1 has been checked
ddmu_dnuda1 = 2.*c.*s./(a^4.*sum_sq.^2) .* (-3.*a0^2.*( -a0 + a2 + a0.*ch0 ...
    -a2.*chin0 ).*(1+cch_sq.*s_sq).*sh + a.*(3.*a.^3.*cth.^2 - (3.*a0.*a1^2 ...
    + a1^3 + a0^2.*(3.*a1+a2) + a0^2.*(a0.*ch0 - a2.*chin0 ) ).*cth.*cch.*...
     (3 + cch_sq.*s_sq) + a^3.*(-1 + cch.^4 .*s_sq ) ) .*sh_sq.*dmu_da1);
                          
% ddmu_dnuda2 has been checked
ddmu_dnuda2 = 2.*c.*s./(a^3.*sum_sq.^2) .* ( - a0^2.*(-1+chin0).*(1+cch_sq.*s_sq).*sh ...
    + (-3.*(3*a0*a1^2 + a1^3 + a0^2*(3*a1+a2) + a0^2.*(a0.*ch0 - a2.*chin0)).*ch ...
    + a^3.*(2+cosh(2.*mu) ) + cch.*(-(3.*a0*a1^2 + a1^3 + a0^2.*(3*a1+a2) ...
    + a0^2*(a0.*ch0 - a2.*chin0)).*cth + a^3*cch ).*s_sq ).*dmu_da2 );

% ddmu_dmu0da1 has been checked
ddmu_dmu0da1 = - a0^3*cch.^3.*sh0.*sum0_sq./(a^4*(1+cch.^2.*s_sq).^2) ...
         .*( 3 + 3.*cch_sq.*s_sq - a/2*(2 + cos(2*nu)-3*cosh(2*mu) ) ...
         .*cth.*cch_sq.*dmu_da1);
     
% ddmu_dmu0da2 has been checked
ddmu_dmu0da2 = -a0^3*cch.*sh0.*sum0_sq.*(cth.*s_sq + 3*ch.*sh).*dmu_da2 ...
                ./(a^3.*sum_sq.^2 );
            
subel.a = a;
subel.nu = nu;
subel.sh = sh;
subel.ch = ch;
subel.sh_sq = sh_sq;
subel.ch_sq = ch_sq;
subel.chin0 = chin0;
subel.sh0 = sh0;
subel.ch0 = ch0;
subel.s = s;
subel.c = c;
subel.c_sq = c_sq;
subel.s_sq = s_sq;
subel.cch_sq = cch_sq;
subel.sum_sq = sum_sq;
subel.sum0_sq = sum0_sq;
subel.root_sum_sq = root_sum_sq;
subel.root_sum0_sq = root_sum0_sq;
subel.dmu_dmu0 = dmu_dmu0;
subel.dmu_dnu  = dmu_dnu ;
subel.dmu_da1  = dmu_da1 ;
subel.dmu_da2 = dmu_da2;
subel.ddmu_dnuda1 = ddmu_dnuda1;
subel.ddmu_dnuda2  = ddmu_dnuda2 ;
subel.ddmu_dmu0da1 = ddmu_dmu0da1;
subel.ddmu_dmu0da2 = ddmu_dmu0da2;
subel.mu = mu;
subel.mu_in = subel.mu(1,:);
subel.dmu_dnu_in = subel.dmu_dnu(1,:);


% Call backs
% a = subel.a;
% nu = subel.nu;
% sh = subel.sh;
% ch = subel.ch;
% sh_sq = subel.sh_sq;
% ch_sq = subel.ch_sq;
% chin0 = subel.chin0;
% sh0 = subel.sh0;
% ch0 = subel.ch0;
% s = subel.s;
% c = subel.c;
% c_sq = subel.c_sq;
% s_sq = subel.s_sq;
% cch_sq = subel.cch_sq;
% sum_sq = subel.sum_sq;
% sum0_sq = subel.sum0_sq;
% root_sum_sq = subel.root_sum_sq;
% root_sum0_sq = subel.root_sum0_sq;
% dmu_dmu0 = subel.dmu_dmu0;
% dmu_dnu  = subel.dmu_dnu ;
% dmu_da1  = subel.dmu_da1 ;
% dmu_da2 = subel.dmu_da2;
% ddmu_dnuda1 = subel.ddmu_dnuda1;
% ddmu_dnuda2  = subel.ddmu_dnuda2 ;
% ddmu_dmu0da1 = subel.ddmu_dmu0da1;
% ddmu_dmu0da2 = subel.ddmu_dmu0da2;
% mu = subel.mu;
% mu_in = subel.mu(1,:);
% dmu_dnu_in = subel.dmu_dnu(1,:);
% muvec0 = subel.muvec0;








end

