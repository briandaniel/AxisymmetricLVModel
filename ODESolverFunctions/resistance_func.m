function [ R_AOV, Rmv ] = resistance_func( Rmv_open, R_AOV_open,  Rmv_close, R_AOV_close, ...
                                smoother, P_A, P_AO, P, R_AO )
%RESISTANCE_FUNC Summary of this function goes here
%   Detailed explanation goes here


Rmv = Rmv_close - ( Rmv_close - Rmv_open) / (1 + exp(- smoother*(P_A - P) ) );

f_aov = @(x) -x + R_AOV_close - ( R_AOV_close - R_AOV_open)/...
            ( 1 + exp(- smoother*(P - ( P*R_AO + P_AO*x )/(R_AO + x) ) ) );
R_AOV = fzero(f_aov,1);

end

