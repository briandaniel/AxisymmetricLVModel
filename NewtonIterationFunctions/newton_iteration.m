function [ f, J, Rmv, Raov ] = newton_iteration( vara, param, X, t)
%NEWTON_SOLVER Summary of this function goes here
%   Detailed explanation goes here

f1 = f1_calc(X, vara, param);
f2 = f2_calc(X, vara, param);
f3 = f3_calc(X, vara, param);
f4 = f4_calc(X, vara, param, t);
f5 = f5_calc(X, vara, param, t);
f = [f1, f2, f3, f4, f5]';

dX = 0.0001;
Xdx4 = [ X(1), X(2), X(3), X(4) + dX, X(5) ];
Xdx5 = [ X(1), X(2), X(3), X(4), X(5) + dX ];

% Calculate derivatives for jacobian
df1_dX1 = vara.alpha_a1;
df1_dX2 = vara.beta_a1;
df1_dX3 = vara.gamma_a1;
df1_dX4 = -vara.eta1;

df2_dX1 = vara.alpha_a2;
df2_dX2 = vara.beta_a2;
df2_dX3 = vara.gamma_a2;
df2_dX4 = -vara.eta2;
        
df3_dX1 = vara.alpha_a3;
df3_dX2 = vara.beta_a3;
df3_dX3 = vara.gamma_a3;
     
df4_dX1 = vara.dV_da1;
df4_dX2 = vara.dV_da2;

df4_dX4 = ( f4_calc(Xdx4, vara, param, t) - f4 )/dX;
df4_dX5 = ( f4_calc(Xdx5, vara, param, t) - f4 )/dX;

A = f4_calc(Xdx5, vara, param, t);
B = f5_calc(Xdx5, vara, param, t);

df5_dX4 = ( f5_calc(Xdx4, vara, param, t) - f5 )/dX;
df5_dX5 = ( f5_calc(Xdx5, vara, param, t) - f5 )/dX;

J = [ df1_dX1,  df1_dX2,  df1_dX3,  df1_dX4,    0;...
      df2_dX1,  df2_dX2,  df2_dX3,  df2_dX4,    0;...
      df3_dX1,  df3_dX2,  df3_dX3,  0,          0;...
      df4_dX1,  df4_dX2,  0,        df4_dX4,    df4_dX5;...
      0,        0,        0,        df5_dX4,    df5_dX5];
  
Raov = Raov_calc(X, vara, param);
Rmv = Rmv_calc(X, vara, param);

end

function f1 = f1_calc(X, vara, param)

    f1 = vara.alpha_a1*X(1) + vara.beta_a1*X(2) + vara.gamma_a1*X(3) ...
            + vara.kappa_a1 - vara.eta1*X(4);

end

function f2 = f2_calc(X, vara, param)

    f2 = vara.alpha_a2*X(1) + vara.beta_a2*X(2) + vara.gamma_a2*X(3) ...
            + vara.kappa_a2 - vara.eta2*X(4);

end

function f3 = f3_calc(X, vara, param)

    f3 = vara.alpha_a3*X(1) + vara.beta_a3*X(2) + vara.gamma_a3*X(3) ...
            + vara.kappa_a3 ;

end

function f4 = f4_calc(X, vara, param, t)
    
    Rmv = Rmv_calc(X, vara, param, t);
    Raov = Raov_calc(X, vara, param);

    f4 = vara.dV_da1*X(1) + vara.dV_da2*X(2) - (vara.Pla - X(4))/Rmv ...
        + (X(4) - X(5))/Raov;

end

function f5 = f5_calc(X,vara,param,t)

    Raov = Raov_calc(X, vara, param);
    f5 = X(5)*(param.Rpao + Raov) - X(4)*param.Rpao - vara.Pao*Raov ;
    
end

function Rmv = Rmv_calc(X, vara, param, t)

    Rmv = param.Rmv_close - ( param.Rmv_close - param.Rmv_open) / ...
                            (1 + exp(- param.smoother*( vara.Pla  - X(4)) ) );

end

function Raov = Raov_calc(X, vara, param)

    Raov = param.Raov_close - ( param.Raov_close - param.Raov_open) / ...
                            (1 + exp(- param.smoother*( X(4)  - X(5)) ) );

end


























