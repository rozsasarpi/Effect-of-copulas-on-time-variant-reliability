function I = innerS(r, e, fRES) 
    I = integral(@(s) fRES(r, e, s), r-e, 6);
end