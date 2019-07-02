function [ phi_a ] = phi_a_calc( t, param )
%PHI_A_CALC Summary of this function goes here
%   Detailed explanation goes here

    phi_a = 0;
    if ( abs(t - param.Tc + param.Tca_shift) <= param.Tca/2 )
        
        phi_a = cos( pi/param.Tca*(t - param.Tc + param.Tca_shift ) ).^2 ;
        
    elseif ( abs(t + param.Tca_shift ) <= param.Tca/2 )
        
        phi_a = cos( pi/param.Tca*(t + param.Tca_shift ) ).^2 ;
        
    end
        

end

