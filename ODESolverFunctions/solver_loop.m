function [ Y, tvec, vara, Xnew, eps_fed] = solver_loop( t0, t1, y0, param, Xguess, eps_fed)
%SOLVER_LOOP Summary of this function goes here
%   Detailed explanation goes here

    dt = param.dt0;
    t = t0;

    tvec = []; 
    vara = [];
    y = y0;
    yn = y; 
    P = [];
       
    k = 1;
    while t <= param.Tc
      

        display(strcat('Computing step t = ', num2str(t)))   

        
        Y(k,:) = y';
        tvec = [tvec;t];
        [ y, Xguess, otherVara, eps_fed] = solver_step( @ode_func, t, yn, dt, param, Xguess, eps_fed);
        yn = y;
        t = t+dt;
        k = k+1;
      
        vara = [vara; otherVara];
        
    end

    y'
    Xnew = Xguess;
    
end












