function [ Sv_ai ] ...
 = viscous_stress( kv, detF, Cinv, dE_dai )

[Nmu,Nnu,~,~] = size(Cinv);
Sv_ai = zeros(size(dE_dai));

for m = 1:3
    for k = 1:3
        Sv_ai(:,:,1,1) = Sv_ai(:,:,1,1) + ...
                        Cinv(:,:,1,m).*dE_dai(:,:,m,k).*Cinv(:,:,k,1);
        Sv_ai(:,:,1,2) = Sv_ai(:,:,1,2) + ...
                        Cinv(:,:,1,m).*dE_dai(:,:,m,k).*Cinv(:,:,k,2); 
        Sv_ai(:,:,1,3) = Sv_ai(:,:,1,3) + ...
                        Cinv(:,:,1,m).*dE_dai(:,:,m,k).*Cinv(:,:,k,3);  

        Sv_ai(:,:,2,1) = Sv_ai(:,:,2,1) + ...
                        Cinv(:,:,2,m).*dE_dai(:,:,m,k).*Cinv(:,:,k,1);
        Sv_ai(:,:,2,2) = Sv_ai(:,:,2,2) + ...
                        Cinv(:,:,2,m).*dE_dai(:,:,m,k).*Cinv(:,:,k,2); 
        Sv_ai(:,:,2,3) = Sv_ai(:,:,2,3) + ...
                        Cinv(:,:,2,m).*dE_dai(:,:,m,k).*Cinv(:,:,k,3);  

        Sv_ai(:,:,3,1) = Sv_ai(:,:,3,1) + ...
                        Cinv(:,:,3,m).*dE_dai(:,:,m,k).*Cinv(:,:,k,1);
        Sv_ai(:,:,3,2) = Sv_ai(:,:,3,2) + ...
                        Cinv(:,:,3,m).*dE_dai(:,:,m,k).*Cinv(:,:,k,2); 
        Sv_ai(:,:,3,3) = Sv_ai(:,:,3,3) + ...
                        Cinv(:,:,3,m).*dE_dai(:,:,m,k).*Cinv(:,:,k,3);  
    end
end



for m = 1:3
    for k = 1:3    
        Sv_ai(:,:,k,m) = Sv_ai(:,:,k,m).*(2*kv).*detF;
    end
end


end


