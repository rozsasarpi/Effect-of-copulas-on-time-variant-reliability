clearvars
% close all
clc

tt      = [0, 2, 5, 10, 15, 20];
xmin    = 100;
xmax    = 13500;
mm      = 300;

xx      = linspace(xmin, xmax, mm);
nn      = length(tt);

for ii = 1:nn
    close all
    t = tt(ii);
    f = nan(mm,1);
    tic
    for jj = 1:mm
        x = xx(jj);
        Fn = @(x) normcdf(-form_wrapper(x, t));
        f(jj) = cfd(Fn, x);       
    end
    toc
    
%     load(['ffAA_',num2str(t),'.mat'])
%     f = f(:);
%     semilogy(xx, f)
%     hold on
%     semilogy(xxAA, ffAA, '--')
    
    xxAA = xx;
    ffAA = f;
    save(['ffBB_0.05_',num2str(t),'.mat'], 'ffAA', 'xxAA')
end


