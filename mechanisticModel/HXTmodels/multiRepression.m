function handle= multiRepression(k1,k2,k3)
   
handle=@assembly;
    
    function final=assembly(r1, r2, r3)
    
       final= (1+denom(r1,k1)+denom(r2,k2)+denom(r3,k3))^-1;
    end
end

        
        
        