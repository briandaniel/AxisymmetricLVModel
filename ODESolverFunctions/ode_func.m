function [dy, X, otherVara, eps_fed] = ode_func(t,Y,param,Xguess, eps_fed)
% Generates left hand side matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load parameters
Nmu = param.Nmu;
Nnu = param.Nnu;
nu_up = param.nu_up;
Tc = param.Tc;
Ta = param.Ta;
a0 = param.a0;
muin0 = param.muin0;
muout0 = param.muout0;
km = param.km;
kav = param.kav;
kp = param.kp;
m = param.m;
kv = param.kv;
cvec = param.cvec;
kd = param.kd;


% load variables
if strcmp( param.odeSystem, 'full')
    Ppa = Y(1);
    Fpa = Y(2);
    Ppp = Y(3);
    Fsa = Y(4);
    Psp = Y(5);
    Pev = Y(6);
    Pra = Y(7);
    Ppv = Y(8);
    Vrv = Y(9);
    Psa = Y(10);
    Vla = Y(11);
    Pao = Y(12);
    a1 = Y(13);
    a2 = Y(14);
    a3 = Y(15);
elseif strcmp( param.odeSystem, 'circOnly')
    Ppa = Y(1);
    Fpa = Y(2);
    Ppp = Y(3);
    Fsa = Y(4);
    Psp = Y(5);
    Pev = Y(6);
    Pra = Y(7);
    Ppv = Y(8);
    Vrv = Y(9);
    Psa = Y(10);
    Pla = Y(11);
    Pao = Y(12);
    a1 = Y(13);
    a2 = Y(14);
    a3 = Y(15);
elseif strcmp( param.odeSystem, 'lvOnly')
    Pla = Y(1);
    Pao = Y(2);
    a1 = Y(3);
    a2 = Y(4);
    a3 = Y(5);
    
    Psa = param.Psa_const;
    Ppv = param.Ppv_const;
else 
    display('Error, no valid ode system chosen')
    return
end
    
dy = zeros(size(Y));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate new mu values
muvec0 = linspace(muin0,muout0,Nmu);
nu_vec = linspace(nu_up,pi,Nnu);

mu = zeros(Nmu,Nnu);
for k=1:Nmu
    mu(k,:) = muCalc( a0, a1, a2, muvec0(k), nu_vec, muin0);
end

subel = subsub_elements( mu, nu_vec, a0, a1, a2, muvec0, muin0 );
subel.muvec0 = muvec0;

% Test deformation tensor
[ F, dF_da1, dF_da2, dF_da3 ] ...
    = deformation_tensor_derivatives( subel, param, a1, a2, a3);

% Test cauchy green defrormation tensor
[ C ] = cauchy_tensor( F );
         
% Test green strain tensor
[ E, dE_da1, dE_da2, dE_da3 ] ...
             = green_strain_derivatives( F, dF_da1, dF_da2, dF_da3, C );
                                   
% Generate rotation matrix
[Q, omega] = initial_rotation_matrix( subel,param );

% generate determinant of F
[ detF ] = determinantF( F );

% Generate Cinv
[ Cinv ] = inverse_cauchy_green_deformation( C );
         
% Generate Sv the viscous stress
[ Sv_a1 ] = viscous_stress( kv, detF, Cinv, dE_da1 );
[ Sv_a2 ] = viscous_stress( kv, detF, Cinv, dE_da2 );
[ Sv_a3 ] = viscous_stress( kv, detF, Cinv, dE_da3 );

% Generate Se the elastic stress
[ Se ] = elastic_stress ( C, E, cvec, Q );


% Generate derivatives of ell
[ ell, ell0, dell_da1, dell_da2, dell_da3 ]  = ...
            ell_derivative( subel, param, a0, a1, a2, a3, omega );

% Compute end-diastolic fiber strain        
if ( t == param.t_init || t == param.Tc )
    display('Computing end diastolic strain')
    eps_fed = 0.5*(ell.^2./ell0.^2 - 1);
end

At = activation_func( t, Ta, Tc, eps_fed, kd );

% Test active fiber stress function  
[ Sf_const ] = active_fiber_stress_constant( At, km, kp, m, ...
                                     ell, ell0, eps_fed, Q, F );
[ Sf_a1 ] = active_fiber_stress( At, kav, kp, m, ...
                                     ell, ell0, eps_fed, dell_da1, Q, F );
[ Sf_a2 ] = active_fiber_stress( At, kav, kp, m, ...
                                     ell, ell0, eps_fed, dell_da2, Q, F );
[ Sf_a3 ] = active_fiber_stress( At, kav, kp, m, ...
                                     ell, ell0, eps_fed, dell_da3, Q, F );

% Traction integral
[eta1, eta2] = traction_integral( subel, param, nu_vec, a0, a1, a2  );

% Left hand equilibrium integral
[ alpha_a1, beta_a1, gamma_a1, kappa_a1 ] ...
    = equilibrium_integral( a0, a1, muvec0, nu_vec, dE_da1, Sv_a1, Sv_a2, Sv_a3, ...
                             Sf_const, Sf_a1, Sf_a2, Sf_a3, Se);
[ alpha_a2, beta_a2, gamma_a2, kappa_a2 ] ...
    = equilibrium_integral( a0, a1, muvec0, nu_vec, dE_da2, Sv_a1, Sv_a2, Sv_a3, ...
    Sf_const, Sf_a1, Sf_a2, Sf_a3, Se);
[ alpha_a3, beta_a3, gamma_a3, kappa_a3 ] ...
    = equilibrium_integral( a0, a1, muvec0, nu_vec, dE_da3, Sv_a1, Sv_a2, Sv_a3, ...
    Sf_const, Sf_a1, Sf_a2, Sf_a3, Se);

% Calculate volume derivatives
[ dV_da1, dV_da2 ] = volume_derivatives( subel,param);

% Calculate LV volume
[ Vlv ] = cavity_volume_calc(subel.mu_in, nu_vec, subel.dmu_dnu_in, subel.a);

% If varying elastance model is used for atrium, calculate the pressure in
% the left atrium
if strcmp( param.odeSystem, 'full')
    % Calculate left atrial pressure
    [ Pla ] = left_atrium_calc( Vla, t, param );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run newton iteration to solve nonlinear system
vara = struct( 'alpha_a1', alpha_a1, 'alpha_a2', alpha_a2, 'alpha_a3', alpha_a3, ... 
               'beta_a1', beta_a1, 'beta_a2', beta_a2, 'beta_a3', beta_a3, ...
               'gamma_a1', gamma_a1, 'gamma_a2', gamma_a2, 'gamma_a3', gamma_a3, ...
               'kappa_a1', kappa_a1, 'kappa_a2', kappa_a2, 'kappa_a3', kappa_a3, ...
               'Ppv', Ppv, 'Psa', Psa, 'Pao', Pao, 'Pla', Pla, ...
               'a1', a1, 'a2', a2, 'a3', a3, ...
               'dV_da1', dV_da1, 'dV_da2', dV_da2, 'eta1', eta1, 'eta2', eta2  );       
[ da1_dt, da2_dt, da3_dt, Plv, Ppao, X, Rmv, Raov] = newton_solver( vara, param, Xguess, t );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Calculate lumped parameter model %%%%%%%%%%%%%%%%%
% Compliances (ml/mmHg)
Cla = param.Cla;
Csa = param.Csa;
Csp = param.Csp;
Cep = param.Cep;
Csv = param.Csv;
Cev = param.Cev;
Cpa = param.Cpa;
Cpp = param.Cpp;
Cpv = param.Cpv;
Cao = param.Cao;

% Zero pressure volumes (ml)
Vu_sa = param.Vu_sa;
Vu_sp = param.Vu_sp;
Vu_ep = param.Vu_ep;
Vu_sv = param.Vu_sv;
Vu_ev = param.Vu_ev;
Vu_pa = param.Vu_pa;
Vu_pp = param.Vu_pp;
Vu_pv = param.Vu_pv;
Vu_la = param.Vu_la;
    
% Hydraulic resistance (mmHg*s/ml)
Rsa = param.Rsa;
Rsp = param.Rsp;
Rep = param.Rep;
Rsv = param.Rsv;
Rev = param.Rev;
Rpa = param.Rpa;
Rpp = param.Rpp;
Rpv = param.Rpv;
Rao = param.Rao;
Rpao = param.Rpao;

% Right heart parameters
Cra = param.Cra ;
Vu_ra = param.Vu_ra;
Rra = param.Rra;
P0_rv = param.P0_rv ;
ke_rv = param.ke_rv;
Vu_rv = param.Vu_rv ;
Emax_rv = param.Emax_rv;
kr_rv = param.kr_rv ;

% Other parameters
Lpa = param.Lpa;
Lsa = param.Lsa;
Vt = param.Vt0;

% Calculate lumped model if circulatory system is being used
if strcmp( param.odeSystem, 'full') 
   
    dy(13) = da1_dt;
    dy(14) = da2_dt;
    dy(15) = da3_dt;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dVla_dt = (Ppv - Pla)/Rpv - (Pla - Plv)/Rmv;         
    dPao_dt = 1/Cao*( (Ppao-Pao)/Rpao - (Pao - Psa)/Rao) ;
    dy(11) = dVla_dt;
    dy(12) = dPao_dt; 
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Left atrial contraction
    [ phi_a ] = phi_a_calc( t, param );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Volume conservation calculations
    Vu = Vu_sa + Vu_sp + Vu_ep + Vu_sv + Vu_ev + Vu_ra + Vu_pa + Vu_pp + ...
         Vu_pv;
    Psv =  1/Csv*(Vt - Csa*Psa - (Csp + Cep)*Psp - Cev*Pev - Cra*Pra ...
            - Vrv - Cpa*Ppa - Cpp*Ppp - Cpv*Ppv - Vla - Vlv - Vu );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    % Right ventricle contraction
    % Activation function calculation
    phi = 0;
    if t < Ta 
        phi = sin(pi*t/Ta);
    end

    Pmax_rv = phi*Emax_rv*(Vrv - Vu_rv) + (1-phi)*P0_rv*(exp(ke_rv*Vrv) - 1);
    Rrv = kr_rv*Pmax_rv;

    % Pulmonary valve
    F_or = 0;
    if (Pmax_rv > Ppa)
        F_or = (Pmax_rv - Ppa)/Rrv;
    end

    Prv = Pmax_rv - Rrv*F_or;

    % tricuspid valve 
    Fir = 0;
    if (Pra > Prv)
        Fir = (Pra - Prv)/Rra;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % d/dt [ P_pa ]
    dy(1) = 1/Cpa*(F_or - Fpa);

    % d/dt [ Fpa ]
    dy(2) = 1/Lpa*(Ppa - Ppp - Rpa*Fpa);

    % d/dt [ Ppp ]
    dy(3) = 1/Cpp*(Fpa - (Ppp-Ppv)/Rpp );

    % d/dt [ Fsa ]
    dy(4) = 1/Lsa*(Psa - Psp - Rsa*Fsa);

    % d/dt [ Psp ]
    dy(5) = 1/(Csp + Cep)*(Fsa - (Psp - Psv)/Rsp - (Psp - Pev)/Rep);

    % d/dt [ Pev ]
    dy(6) = 1/Cev*( (Psp - Pev)/Rep - (Pev - Pra)/Rev );

    % d/dt [ Pra ]
    dy(7) = 1/Cra * ( (Psv - Pra)/Rsv + (Pev - Pra)/Rev - Fir);

    % d/dt [ Ppv ]
    dy(8) = 1/Cpv*( (Ppp - Ppv)/Rpp - (Ppv - Pla)/Rpv);

    % d/dt [ Vrv ]
    dy(9) = Fir - F_or;

    % d/dt [ Psa ]
    dy(10) = 1/Csa*( (Pao - Psa)/Rao - Fsa );

    otherVara(1) = Plv;
    otherVara(2) = Rmv;
    otherVara(3) = Raov;
    otherVara(4) = Prv;
    otherVara(5) = Vlv;
    otherVara(6) = Vrv;
    otherVara(7) = Pla;
    otherVara(8) = phi;
    otherVara(9) = phi_a;

elseif strcmp( param.odeSystem, 'circOnly')
   
    dy(13) = da1_dt;
    dy(14) = da2_dt;
    dy(15) = da3_dt;
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dPla_dt =  1/param.Cla*((Ppv - Pla)/Rpv - (Pla - Plv)/Rmv);
    dy(11) = dPla_dt;
    
    dPao_dt = 1/Cao*( (Ppao-Pao)/Rpao - (Pao - Psa)/Rao) ;
    dy(12) = dPao_dt; 
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Left atrial contraction
    [ phi_a ] = phi_a_calc( t, param );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Volume conservation calculations
    Vu = Vu_sa + Vu_sp + Vu_ep + Vu_sv + Vu_ev + Vu_ra + Vu_pa + Vu_pp + ...
     Vu_pv + Vu_la;
    Psv =  1/Csv*(Vt - Csa*Psa - (Csp + Cep)*Psp - Cev*Pev - Cra*Pra ...
        - Vrv - Cpa*Ppa - Cpp*Ppp - Cpv*Ppv - Cla*Pla - Vlv - Vu );


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    % Right ventricle contraction
    % Activation function calculation
    phi = 0;
    if t < Ta 
        phi = sin(pi*t/Ta);
    end

    Pmax_rv = phi*Emax_rv*(Vrv - Vu_rv) + (1-phi)*P0_rv*(exp(ke_rv*Vrv) - 1);
    Rrv = kr_rv*Pmax_rv;

    % Pulmonary valve
    F_or = 0;
    if (Pmax_rv > Ppa)
        F_or = (Pmax_rv - Ppa)/Rrv;
    end

    Prv = Pmax_rv - Rrv*F_or;

    % tricuspid valve 
    Fir = 0;
    if (Pra > Prv)
        Fir = (Pra - Prv)/Rra;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % d/dt [ P_pa ]
    dy(1) = 1/Cpa*(F_or - Fpa);

    % d/dt [ Fpa ]
    dy(2) = 1/Lpa*(Ppa - Ppp - Rpa*Fpa);

    % d/dt [ Ppp ]
    dy(3) = 1/Cpp*(Fpa - (Ppp-Ppv)/Rpp );

    % d/dt [ Fsa ]
    dy(4) = 1/Lsa*(Psa - Psp - Rsa*Fsa);

    % d/dt [ Psp ]
    dy(5) = 1/(Csp + Cep)*(Fsa - (Psp - Psv)/Rsp - (Psp - Pev)/Rep);

    % d/dt [ Pev ]
    dy(6) = 1/Cev*( (Psp - Pev)/Rep - (Pev - Pra)/Rev );

    % d/dt [ Pra ]
    dy(7) = 1/Cra * ( (Psv - Pra)/Rsv + (Pev - Pra)/Rev - Fir);

    % d/dt [ Ppv ]
    dy(8) = 1/Cpv*( (Ppp - Ppv)/Rpp - (Ppv - Pla)/Rpv);

    % d/dt [ Vrv ]
    dy(9) = Fir - F_or;

    % d/dt [ Psa ]
    dy(10) = 1/Csa*( (Pao - Psa)/Rao - Fsa );

    otherVara(1) = Plv;
    otherVara(2) = Rmv;
    otherVara(3) = Raov;
    otherVara(4) = Prv;
    otherVara(5) = Vlv;
    otherVara(6) = Vrv;
    otherVara(7) = Pla;
    otherVara(8) = phi;
    otherVara(9) = phi_a;  
       
elseif strcmp( param.odeSystem, 'lvOnly')

    dPla_dt =  1/param.Cla*((Ppv - Pla)/Rpv - (Pla - Plv)/Rmv);
    dy(1) = dPla_dt;
    
    dPao_dt = 1/Cao*( (Ppao-Pao)/Rpao - (Pao - Psa)/Rao) ;
    dy(2) = dPao_dt; 
       
    dy(3) = da1_dt;
    dy(4) = da2_dt;
    dy(5) = da3_dt;
    
    otherVara(1) = Plv;
    otherVara(2) = Rmv;
    otherVara(3) = Raov;
    otherVara(4) = Vlv;
    otherVara(5) = mean(mean(At));
else 
    display('Error, no valid ode system chosen')
    return
end


end























