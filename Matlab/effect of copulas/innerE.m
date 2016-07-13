function I = innerE(r, fRE) 
    I = integral(@(e) fRE(r,e), r, 6);
end