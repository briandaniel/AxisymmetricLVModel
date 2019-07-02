clear all
close all
clc



% Load parameters and set initial values
run parameter_values.m

if strcmp( odeSystem, 'full')
    Y0 = [Ppa_init, Fpa_init, Ppp_init, Fsa_init, Psp_init, Pev_init, ...
        Pra_init, Ppv_init, Vrv_init, Psa_init, Vla_init, Pao_init, ...
        a1_init, a2_init, a3_init]';
elseif strcmp( odeSystem, 'circOnly')
    Y0 = [Ppa_init, Fpa_init, Ppp_init, Fsa_init, Psp_init, Pev_init, ...
        Pra_init, Ppv_init, Vrv_init, Psa_init, Pla_init, Pao_init, ...
        a1_init, a2_init, a3_init]';
elseif strcmp( odeSystem, 'lvOnly')
    Y0 = [Pla_init, Pao_init, a1_init, a2_init, a3_init]';  
else 
    display('Error, no valid ode system chosen')
    return
end


% Run preloop
tic

% Calculate initial epsilon_fed (typically 0)
eps_fed = 0;

% Run preloop before atrial contraction starts if necessary
% if strcmp( odeSystem, 'full')
%     [ y, t, otherVara, Xnew, eps_fed] =  solver_loop( Tc - Tca_shift - Tca, Tc, Y0, param, Xguess, eps_fed);
%     Y0 = y(end,:)';
%     T = t - t(1);
%     Y = y;
%     otherVara_vec = otherVara;
%     tend = T(end);        
% else
    t = t_init;
    T = [];
    Y = [];
    otherVara_vec = [];
    tend = 0;
% end
   
for N = 1:Ncycles
    
    display(strcat('Currently computing cycle: ', num2str(N) ) )

    [ y, t, otherVara, Xnew, eps_fed] =  solver_loop( t_init, Tc, Y0, param, Xguess, eps_fed);
    Y0 = y(end,:)';
    T = [T; (t + (N-1)*Tc + tend)];
    Y = [Y; y];
    otherVara_vec = [otherVara_vec; otherVara];
    
    % load variables
    if strcmp( param.odeSystem, 'full') || strcmp( param.odeSystem, 'circOnly')
        a1 = y(13);
        a2 = y(14);
        a3 = y(15);
    elseif strcmp( param.odeSystem, 'lvOnly')
        a1 = y(3);
        a2 = y(4);
        a3 = y(5);
    else 
        display('Error, no valid ode system chosen')
        return
    end

end
   
toc 

simTime = toc;



%% Visualize data

run visualize_v20.m
























