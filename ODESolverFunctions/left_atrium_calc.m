function [ Pla ] = left_atrium_calc( Vla, t, param )
%LEFT_ATRIUM_CALC Summary of this function goes here
%   Detailed explanation goes here

    [ phi_a ] = phi_a_calc( t, param );

%     ya = -12*phi_a + 0.5;
    ya = phi_a*20;
    
    Ela = (param.Ela_max - param.Ela_min)*ya + param.Ela_min;
    Vla_rest = (1-ya)*(param.Vla_rd - param.Vla_rs) + param.Vla_rs;

    Pla = Ela*(Vla - Vla_rest);
 
end

