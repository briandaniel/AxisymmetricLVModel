function [ B ] = normalMatrix( A )
%NORMALMATRIX Summary of this function goes here
%   Detailed explanation goes here

[Nmu, Nnu, m, n] = size(A);

B = zeros(m,n,Nmu,Nnu);

for i = 1:Nmu
    for j = 1:Nnu
        for k = 1:m
            for l = 1:n
                
                B(k,l,i,j) = A(i,j,k,l);
            end
        end
    end
end
        
end

