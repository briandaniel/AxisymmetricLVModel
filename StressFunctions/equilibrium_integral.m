function [ alpha_ai, beta_ai, gamma_ai, kappa_ai ] ...
        = equilibrium_integral( a0, a1, mu0_vec, nu_vec, dE_dai, ...
                Sv_a1, Sv_a2, Sv_a3, Sf_const, Sf_a1, Sf_a2, Sf_a3, Se)
%EQUILIBRIUM_INTEGRAL Summary of this function goes here
%   Detailed explanation goes here

Nmu = length(mu0_vec);
Nnu = length(nu_vec);
nu =  repmat( nu_vec, Nmu, 1);
mu0 = repmat( mu0_vec', 1, Nnu );

gmu0 = a0*sqrt(sinh(mu0).^2 + sin(nu).^2);
gnu0 = gmu0;
gphi0 = a0*sinh(mu0).*sin(nu);

integratingFactor = gmu0.*gnu0.*gphi0;

% Calculate constant term coefficient: kappa
kappa_integrand = zeros(Nmu,Nnu);
for k = 1:3
    for j = 1:3
        kappa_integrand = kappa_integrand + dE_dai(:,:,k,j).*(Se(:,:,k,j) + Sf_const(:,:,k,j));
    end
end
kappa_integrand = kappa_integrand.*integratingFactor;
kappa_integral = doubleIntegral( mu0_vec, nu_vec, kappa_integrand );
kappa_ai = 2*pi*kappa_integral;


% Calculate a_1(dot) term coefficient: alpha
alpha_integrand = zeros(Nmu,Nnu);
for k = 1:3
    for j = 1:3
        alpha_integrand = alpha_integrand + dE_dai(:,:,k,j).*(Sv_a1(:,:,k,j) + Sf_a1(:,:,k,j));
    end
end
alpha_integrand = alpha_integrand.*integratingFactor;
alpha_integral = doubleIntegral(  mu0_vec, nu_vec, alpha_integrand );
alpha_ai = 2*pi*alpha_integral;


% Calculate a_2(dot) term coefficient: beta
beta_integrand = zeros(Nmu,Nnu);
for k = 1:3
    for j = 1:3
        beta_integrand = beta_integrand + dE_dai(:,:,k,j).*(Sv_a2(:,:,k,j) + Sf_a2(:,:,k,j));
    end
end
beta_integrand = beta_integrand.*integratingFactor;
beta_integral = doubleIntegral( mu0_vec, nu_vec, beta_integrand );
beta_ai = 2*pi*beta_integral;

% Calculate a_3(dot) term coefficient: gamma
gamma_integrand = zeros(Nmu,Nnu);
for k = 1:3
    for j = 1:3
        gamma_integrand = gamma_integrand + dE_dai(:,:,k,j).*(Sv_a3(:,:,k,j) + Sf_a3(:,:,k,j));
    end
end
gamma_integrand = gamma_integrand.*integratingFactor;
gamma_integral = doubleIntegral( mu0_vec, nu_vec, gamma_integrand );
gamma_ai = 2*pi*gamma_integral;




end




