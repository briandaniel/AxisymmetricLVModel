function [ F, dF_da1, dF_da2, dF_da3 ] ...
 = deformation_tensor_derivatives( subel, param, a1, a2, a3);

% deformation_tensor_derivatives
% inputs
%  1. mu_vec should be a row vector
%  2. nu_vec should be a row vector
%
% Calculates the a_i derivatives of the deformation gradient tensor
% The matrices are arranged with
%       changing mu values across the rows
%       changing nu values down the columns


% Call backs
a = subel.a;
sh = subel.sh;
ch = subel.ch;
sh_sq = subel.sh_sq;
sh0 = subel.sh0;
s_sq = subel.s_sq;
root_sum_sq = subel.root_sum_sq;
root_sum0_sq = subel.root_sum0_sq;
dmu_dmu0 = subel.dmu_dmu0;
dmu_dnu  = subel.dmu_dnu ;
dmu_da1  = subel.dmu_da1 ;
dmu_da2 = subel.dmu_da2;
ddmu_dnuda1 = subel.ddmu_dnuda1;
ddmu_dnuda2  = subel.ddmu_dnuda2 ;
ddmu_dmu0da1 = subel.ddmu_dmu0da1;
ddmu_dmu0da2 = subel.ddmu_dmu0da2;
Nmu = param.Nmu;
Nnu = param.Nnu;
a0 = param.a0;

% Generate matrices
F = zeros(Nmu, Nnu, 3, 3);
dF_da1 = zeros(Nmu, Nnu, 3, 3);
dF_da2 = zeros(Nmu, Nnu, 3, 3);
dF_da3 = zeros(Nmu, Nnu, 3, 3);


% Calculate subelements
g = a*root_sum_sq./(a0*root_sum0_sq);
dg_da1 = (s_sq + sh_sq + a*ch.*sh.*dmu_da1)...
                    ./(a0*root_sum0_sq.*root_sum_sq);
dg_da2 = (a * ch .* sh .* dmu_da2 )...
                    ./(a0*root_sum0_sq.*root_sum_sq);

                
% F11 calculations
F11 = g.*dmu_dmu0;

dF11_da1 = dg_da1.*dmu_dmu0 + g.*ddmu_dmu0da1;
dF11_da2 = dg_da2.*dmu_dmu0 + g.*ddmu_dmu0da2;
dF11_da3 = 0;


% F12 calculations
F12 = g.*dmu_dnu;

dF12_da1 = dg_da1.*dmu_dnu + g.*ddmu_dnuda1;
dF12_da2 = dg_da2.*dmu_dnu + g.*ddmu_dnuda2;
dF12_da3 = 0;


% F22 calculations
F22 = g;

dF22_da1 = dg_da1;
dF22_da2 = dg_da2;
dF22_da3 = 0;


% F32 calculations
F32 = -a3*a*sh.*s_sq./(a0.*root_sum0_sq);

dF32_da1 = - a3*s_sq.*(sh + a*ch.*dmu_da1)./(a0.*root_sum0_sq);
dF32_da2 = - a*a3*ch.*s_sq.*dmu_da2./(a0.*root_sum0_sq);
dF32_da3 = -a*sh.*s_sq./(a0.*root_sum0_sq);


% F33 calculations
F33 = a*sh./(a0*sh0);

dF33_da1 = (sh + a*ch.*dmu_da1)./(a0.*sh0);
dF33_da2 = a.*ch.*dmu_da2./(a0*sh0);
dF33_da3 = 0;

F(:,:,1,1) = F11;
F(:,:,1,2) = F12;
F(:,:,2,2) = F22;
F(:,:,3,2) = F32;
F(:,:,3,3) = F33;

dF_da1(:,:,1,1) = dF11_da1;
dF_da1(:,:,1,2) = dF12_da1;
dF_da1(:,:,2,2) = dF22_da1;
dF_da1(:,:,3,2) = dF32_da1;
dF_da1(:,:,3,3) = dF33_da1;

dF_da2(:,:,1,1) = dF11_da2;
dF_da2(:,:,1,2) = dF12_da2;
dF_da2(:,:,2,2) = dF22_da2;
dF_da2(:,:,3,2) = dF32_da2;
dF_da2(:,:,3,3) = dF33_da2;

dF_da3(:,:,1,1) = dF11_da3;
dF_da3(:,:,1,2) = dF12_da3;
dF_da3(:,:,2,2) = dF22_da3;
dF_da3(:,:,3,2) = dF32_da3;
dF_da3(:,:,3,3) = dF33_da3;

end


