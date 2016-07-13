function I = innerES(r, fRES) 
    I = integral(@(e) innerS(r, e, fRES), -6, 6, 'ArrayValued', true);
end