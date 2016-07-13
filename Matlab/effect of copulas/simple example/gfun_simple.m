function g = gfun_simple(S1, S2, R, als)
       
    switch als
        case 1
            S = S1;
        case 2
            S = S2;
    end
    
    g  = R -S;
    
end


