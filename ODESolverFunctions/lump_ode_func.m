function [ dy ] = ode_func( t, y, param)
%ODE_FUNC Summary of this function goes here
%   Detailed explanation goes here

display(t)
dy = zeros(length(y),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load parameters
% Compliances (ml/mmHg)
Csa = param.Csa;
Csp = param.Csp;
Cep = param.Cep;
Csv = param.Csv;
Cev = param.Cev;
Cpa = param.Cpa;
Cpp = param.Cpp;
Cpv = param.Cpv;

% Zero pressure volumes (ml)
Vu_sa = param.Vu_sa;
Vu_sp = param.Vu_sp;
Vu_ep = param.Vu_ep;
Vu_sv = param.Vu_sv;
Vu_ev = param.Vu_ev;
Vu_pa = param.Vu_pa;
Vu_pp = param.Vu_pp;
Vu_pv = param.Vu_pv;

% Hydraulic resistance (mmHg*s/ml)
Rsa = param.Rsa;
Rsp = param.Rsp;
Rep = param.Rep;
Rsv = param.Rsv;
Rev = param.Rev;
Rpa = param.Rpa;
Rpp = param.Rpp;
Rpv = param.Rpv;

% Left heart parameters
Cla = param.Cla;
Vu_la = param.Vu_la;
Rla = param.Rla ;
P0_lv = param.P0_lv;
ke_lv = param.ke_lv ;
Vu_lv = param.Vu_lv ;
Emax_lv = param.Emax_lv;
kr_lv = param.kr_lv;

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
ksys = param.ksys;
Tsys_0 = param.Tsys_0;
T = param.T0;
Vt = param.Vt0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dVu_ev_dt = 0; % for simplicity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load y values into variable names for readable coding
Ppa = y(1);
Fpa = y(2);
Ppp = y(3);
Ppv = y(4);
Psa = y(5);
Fsa = y(6);
Psp = y(7);
Pev = y(8);
Pla = y(9);
Vlv = y(10);
xi = y(11);
Pra = y(12);
Vrv = y(13);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% d/dt [ Fpa ]
dy(2) = 1/Lpa*(Ppa - Ppp - Rpa*Fpa);

% d/dt [ Ppp ]
dy(3) = 1/Cpp*(Fpa - (Ppp-Ppv)/Rpp );

% d/dt [ Ppv ]
dy(4) = 1/Cpv*( (Ppp - Ppv)/Rpp - (Ppv - Pla)/Rpv);

% d/dt [ Fsa ]
dy(6) = 1/Lsa*(Psa - Psp - Rsa*Fsa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General calculations
Vu = Vu_sa + Vu_sp + Vu_ep + Vu_sv + Vu_ev + Vu_ra + Vu_pa + Vu_pp + ...
     Vu_pv + Vu_la;
Psv =  1/Csv*(Vt - Csa*Psa - (Csp + Cep)*Psp - Cev*Pev - Cra*Pra ...
        - Vrv - Cpa*Ppa - Cpp*Ppp - Cpv*Ppv - Cla*Pla - Vlv - Vu );
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

% d/dt [ Psp ]
dy(7) = 1/(Csp + Cep)*(Fsa - (Psp - Psv)/Rsp - (Psp - Pev)/Rep);

% d/dt [ Pev ]
dy(8) = 1/Cev*( (Psp - Pev)/Rep - (Pev - Pra)/Rev - dVu_ev_dt );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = xi - floor(xi); % Determines phase of heart cycle
Tsys = Tsys_0 - ksys*1/T; % Length of systole

% Activation function calculation
phi = 0;
if ( u <= Tsys/T )
    phi = (sin(pi*T/Tsys*u))^2;
end

Pmax_lv = phi*Emax_lv*(Vlv - Vu_lv) + (1-phi)*P0_lv*(exp(ke_lv*Vlv) - 1);

Rlv = kr_lv * Pmax_lv

% aortic valve
Fol = 0;
if (Pmax_lv > Psa)
    Fol = (Pmax_lv - Psa)/Rlv;
end

Plv = Pmax_lv - Rlv*Fol;

% atrioventricular (mitral) valve
Fil = 0;
if (Pla > Plv)
    Fil = (Pla - Plv)/Rla;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% d/dt [ Psa ]
dy(5) = 1/Csa*( Fol - Fsa );

% d/dt [ Pla ]
dy(9) = 1/Cla*( (Ppv-Pla)/Rpv - Fil );

% d/dt [ Vlv ]
dy(10) = Fil - Fol;

% d/dt [ xi ]
dy(11) = 1/T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% d/dt [ Pra ]
dy(12) = 1/Cra * ( (Psv - Pra)/Rsv + (Pev - Pra)/Rev - Fir);

% d/dt [ Vrv ]
dy(13) = Fir - F_or;



end

