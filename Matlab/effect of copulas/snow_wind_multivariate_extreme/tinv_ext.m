% tinv() extension for 2 degrees of freedom!
% to evaluate points of the t distribution which gives NaN in matlab built in
% not effective but at least it is working

%=============================
% to evaluate points of the t distribution which gives NaN in matlab built in
% PP = logspace(-50, -15, 2e3);
% % PP = 1-PP;
% 
% For df = 2, there is analytical solution! [http://en.wikipedia.org/wiki/Quantile_function#The_Student.27s_t-distribution]
%
% x = (2*p-1)./sqrt(2*p.*(1-p))
%
% tt = zeros(numel(PP),1);
% for i = 1:numel(PP)
%     tt(i) = fzero(@(x) tcdf(x,2) - PP(i), 0);
% end
% plot(PP, tt, 'o-')
% 
% PP_low = PP;
% tt_low = tt;
% save('tinv_ext_points_low.mat', 'PP_low', 'tt_low')
%=============================

function t = tinv_ext(P, nu)
    load('tinv_ext_points_low.mat', 'PP_low', 'tt_low')
%     load('tinv_ext_points_up.mat', 'PP_up', 'tt_up')
    t           = zeros(size(P));
%     keyboard
    idx_lowout  = P <= 1e-40;
    idx_low     = P < 1e-15 & P > 1e-40;
    idx_up      = P > 1-1e-15;
    idx_mid     = not(idx_low | idx_up | idx_lowout);
    
    t(idx_lowout)= min(tt_low); % WARNING!!
    t(idx_low)   = interp1(PP_low, tt_low, P(idx_low));
    t(idx_mid)   = tinv(P(idx_mid), nu);
%     t(idx_up)    = interp1(PP_up, tt_up, P(idx_up));
    t(idx_up)    = tinv(P(idx_up), nu);
    
%     t = reshape(t, size(P,1), size(P,2));
end
