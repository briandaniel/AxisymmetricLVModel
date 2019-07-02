function [ detF ] = determinantF( F )
%DETERMINANTF Summary of this function goes here
%   Detailed explanation goes here

    
    detF = F(:,:,1,1).*F(:,:,2,2).*F(:,:,3,3);

end

