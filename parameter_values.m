% parameter_values.m
% Assigns values to all relevant model parameters

% Add function folders to the path
addpath('FigureFunctions', 'KinematicFunctions', 'StrainFunctions',...
        'StressFunctions', 'UtilityFunctions','ODESolverFunctions',...
        'Figures', 'FunctionTests', 'FigureFunctions\figsaver', 'NewtonIterationFunctions')
        
    
% Choose which ode system to run    
% odeSystem = 'full';     % Run circulation with left atrium contraction
% odeSystem = 'circOnly'; % Run circulation without left atrium contraction
odeSystem = 'full'; % No circulation or left atrium contraction
    
    
% Physical domain parameters
Nmu = 10;%matched                   % Number of points in the \mu direction
Nnu = 10; %matched                   % Number of points in the \nu direction
nu_up = 1.57;%matched             % Smallest \nu value ( \nu ranges in [ nu_up, pi] )
nu_vec = linspace(nu_up,pi,Nnu);    % \nu vector

% Time (/cycle) parameters
t_init = 0;                         % Initial time value (s)
Ncycles = 5;                        % Number of cardiac cycles
dt0 = 0.002;
dt_cutoff = .02;

newton_kmax = 50;
Xguess = [1, 1, 1, 1, 1]';      % Initial guess for newton iteration

% Ventricular contraction
Tc = 1;                             % Length of full cardiac cycle (s) 
Ta = .45;                            % Length of active contraction (s)

% Atrial contraction
Tca = .23; % Length of atrial contraction
Tca_shift = 0.00;

% Initial structure parameters
a0 = 6.2;  %probably matched (called a_focus?)  % Initial spheroid length (cm)
muin0 = .425; %matched              % Inner wall starting \mu value
muout0 = .645; %matched              % Outer wall starting \mu value

% Muscle fiber parameterization angles
psi_in_b0 = -85 * (pi / 180); %matched   % Angle of fiber wrapping on endocardial surface
psi_out_b0 = 65 * (pi / 180); %matched   % Angle of fiber wrapping on epidcardial surface

% Active stress parameters
km = 50;                            % Initial strength of active fiber stress
kav = .5;                           % Coefficient of the viscous component of the active fiber stress
kp = 5;                             % Coefficient of the end-diastolic response (pre-load dependent active stress)
kd = 10;
m = 0.458;                          % Slope of the stress/strain relationship in the active fibers

% Passive stress parameters
kv = 0.25;                          % Viscous stress coefficient
cvec = [5,.7,.7,.7];                % Elastic stress coefficients
                                    % c1: overall energy multiplier (outside exponent)
                                    % c2: 1/4(bff-bxx), (modifies, principally, the fiber strain)
                                    % c3: 1/4bxx, Modifies the principle strain invariants (orthogonal directions)
                                    % c4: 1/2(bfx-bxx), modifies the cross terms

% Variable resistances of the valves
smoother = 30;                      % Smoothing parameter for open and closing valves (beta in the text)
Rmv_open = .001;                    % Resistance of the open mitral valve
Raov_open = .001;                   % Resistance of the open aortic valve
Rmv_close = 50;                   % Resitance of the closed mitral valve
Raov_close = 50;                  % Resistance of the open mitral valve
Raov_maxerror = Rmv_open*10^(-3);
Raov_guess = Rmv_open;

% Initial values
Ppa_init = 1.7867;
Fpa_init = 6.5;
Ppp_init = 1.7667;
Fsa_init = 5.9;
Psp_init = 11.96;
Pev_init = 0.6533;
Pra_init = 0.6666;
Ppv_init = 1.1067;
Vrv_init = 132.95;
Psa_init = 12;
Vla_init = 41     ;  
Pao_init = 11.96  ;                       
a1_init = 0   ;                 
a2_init = 0   ;                 
a3_init = 0 ;

% Values for substitute if system does not include atrial contraction
Pla_init = 1;
Cla = 60;
Vu_la = 25;

% Static pressures if there is no circulation
Ppv_const = 1;
Psa_const = 8;

%%%%%%%%%%% LUMP PARAMETER VALUES %%%%%%%%%%%%%%%%%%%%
% Compliances of the blood flow model
Cao = 1.35;
Rpao = 0.0001      ;        
Rao = 0.0026667  ;

% Compliances (converted from ml/mmHg)
Cpv = 190.28;
Cra = 234.38;
Csa = 2.11;
Csp = 15.375;
Cep = 12.525;
Csv = 458.32;
Cev = 375;
Cpa = 5.7;
Cpp = 43.5;

% Zero pressure volumes (ml)
Vu_sa =    0.0  ;
Vu_sp =  274.4  ;
Vu_ep =  336.6  ;
Vu_sv = 1121.0  ;
Vu_ev = 1375.0  ;
Vu_pa =    0.0  ;
Vu_pp =  123.0  ;
Vu_pv =  120.0  ;
Vu_ra =   25.0  ;
Vu_rv =   40.8  ;


% Hydraulic resistance (converted from  mmHg*s/ml)
Rsa = 0.00533;
Rsp = 0.44093;
Rep = 0.1876;
Rsv = 0.00506;
Rev = 0.00213;
Rpa = 0.00306;
Rpp = 0.01192;
Rpv = 0.00448;
Rra = 0.000333;
Rla = 0.000333;

% Innertance ( converted from mmHg*ml/s^2)
Lpa = 2.9333e-05 ;
Lsa = 2.4e-05;

% Right heart parameters
P0_rv = 0.2			;
ke_rv = 0.011		;		
Emax_rv =  0.2333   ;
kr_rv = 0.0014		;
ksys = .075;
Vt0 = 5300;

% Atrial contraction parameters
Vla_rd = 24;
Vla_rs = 23;
Ela_max = 0.0782;
Ela_min = 0.0711;


param = struct( 'Nmu', Nmu, 'Nnu', Nnu, 'nu_up', nu_up, 'nu_vec', nu_vec, ...
    't_init', t_init, 'Tc', Tc, 'Ta', Ta, 'a0', a0,...
    'muin0', muin0, 'muout0', muout0, 'psi_in_b0', psi_in_b0,  ...
    'psi_out_b0', psi_out_b0, 'km', km, 'kav', kav, 'kp', kp, 'kd',kd, 'm', m,...
    'kv', kv, 'cvec', cvec, 'Cao', Cao, 'Rpao', Rpao, 'Rao', Rao, ...
    'smoother', smoother, 'Rmv_open', Rmv_open, 'Raov_open', Raov_open,...
    'Rmv_close', Rmv_close, 'Raov_close', Raov_close, 'dt0', dt0, ...
    'dt_cutoff', dt_cutoff, 'Csa',Csa,'Csp',Csp,'Cep',Cep,'Csv',Csv,'Cev',Cev,'Cpa',Cpa,'Cpp',Cpp,'Cpv',Cpv, ...
    'Vu_sa',Vu_sa,'Vu_sp',Vu_sp,'Vu_ep',Vu_ep,'Vu_sv',Vu_sv,'Vu_ev',Vu_ev, ...
    'Vu_pa',Vu_pa,'Vu_pp',Vu_pp,'Vu_pv',Vu_pv,...
    'Rsa',Rsa,'Rsp',Rsp,'Rep',Rep,'Rsv',Rsv,'Rev',Rev,'Rpa',Rpa,'Rpp',Rpp,'Rpv',Rpv, ...
    'Cra',Cra,'Vu_ra',Vu_ra,'Rra',Rra,'P0_rv',P0_rv,'ke_rv',ke_rv,'Vu_rv',Vu_rv,'Emax_rv',Emax_rv,...
    'kr_rv',kr_rv,'Lpa',Lpa,'Lsa',Lsa,'ksys',ksys,'Vt0',Vt0,'Tca_shift', Tca_shift, 'Tca', Tca, ...
    'Vla_rd', Vla_rd, 'Vla_rs', Vla_rs, 'Ela_max', Ela_max, 'Ela_min', Ela_min,...
    'newton_kmax', newton_kmax, 'Raov_maxerror', Raov_maxerror, ...
    'Raov_guess', Raov_guess, 'a1_init', a1_init, 'a2_init', a2_init, 'a3_init', a3_init,...
    'odeSystem', odeSystem, 'Pla_init', Pla_init, 'Cla', Cla, 'Vu_la', Vu_la, ...
    'Psa_const', Psa_const, 'Ppv_const', Ppv_const );


      

































