function [ da1_dt, da2_dt, da3_dt, Plv, Ppao, X, Rmv, Raov] = newton_solver( vara, param, Xguess, t)
%NEWTON_SOLVER Summary of this function goes here
%   Detailed explanation goes here

    X = Xguess;    
    
    Xold = ones(5,1)*10^10;
    fold = Xold;
    f = zeros(5,1);
    k = 0;
    
   
    while max(abs(Xold-X)) > 1e-6 && max(abs(fold-f)) > 1e-10 && k < param.newton_kmax;
        Xold = X;
        fold = f;
        [ f, J, Rmv, Raov ] = newton_iteration( vara, param, X, t);
        X = X - J\f;
        k = k+1;
        
    end

    
    if k == param.newton_kmax
        display('Newton solver did not converge')
         abs(sum(Xold-X))
         
         vara
         
         param
         
         X
    end

    da1_dt = X(1);
    da2_dt = X(2); 
    da3_dt = X(3); 
    Plv = X(4);
    Ppao = X(5);

end

