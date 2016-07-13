function Dn = debye(x, n)
    % Debye function.
    %   Dn = DEBYE(X, n) returns the nth order Debye function, evaluated at X.
    %   X is a positive scalar
    %
    %      (n/x^n) * integral from 0 to x of (t^n/(exp(t)-1)) dt
    %
    
    if abs(x) >= realmin
        Dn = n/x^n*integral(@(t) (t.^n)./(exp(t)-1), 0, x);
    else
        Dn = 1;
    end
    
end