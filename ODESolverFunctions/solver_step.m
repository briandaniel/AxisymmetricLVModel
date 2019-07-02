function [ y, newX, otherVara, eps_fed] = solver_step( ode_func, t, yn, dt, param, Xguess, eps_fed)
%RK4 Summary of this function goes here
%   Detailed explanation goes here

%     [fn, newX] = ode_func( t, yn, param, Xguess );
%     y = yn + dt*(23/12*fn - 4/3*fn_m1 + 5/12*fn_m2);
%     y = yn + dt*( 3/2*fn - 1/2*fn_m1 );
%     y = yn + dt*( fn );


% % Fourth order runge-kutta
%     [k1, newX, otherVara, eps_fed] = ode_func( t, yn, param, Xguess, eps_fed);
% 
%     [k2, newX, otherVara, eps_fed] = ode_func( t + dt/2, yn + dt/2*k1, param, newX, eps_fed);
% 
%     [k3, newX, otherVara, eps_fed] = ode_func( t + dt/2, yn + dt/2*k2, param, newX, eps_fed);
% 
%     [k4, newX, otherVara, eps_fed] = ode_func( t + dt, yn + dt*k3, param, newX, eps_fed);
% 
%     y = yn + dt/6*(k1 + 2*k2 + 2*k3 + k4);
%     
% %     k1
%     k2
%     k3
%     k4
%     

% y
%     pause(1)

%     % Second order runge-kutta
%     [k1, newX, otherVara] = ode_func( t, yn, param, Xguess );
%     [k2, newX, otherVara] = ode_func( t + dt/2, yn + dt/2*k1, param, newX );
% 
%     y = yn + dt*k2;


[dy, newX, otherVara, eps_fed] = ode_func( t, yn, param, Xguess, eps_fed);

y = yn + dt*( dy );
end

