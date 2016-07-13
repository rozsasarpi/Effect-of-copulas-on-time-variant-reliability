% FORM analysis using FERUM
% units: [N], [m]

function [beta, formresults, probdata] = form_wrapper(r, t)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'PROBDATA'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Names of random variables. Default names are 'x1', 'x2', ..., if not explicitely defined.

probdata.name = { 'fy';
                  'b0';
                  'h0';
                  'r';
                  't'};


% Marginal distributions for each random variable
% probdata.marg =  [ (type) (mean) (stdv) (startpoint) (p1) (p2) (p3) (p4) (input_type); ... ];

% s_g's distribution is defined with distribution parameters, input_type = 1!
probdata.marg =  [  2   240    0.1*240   240    NaN      NaN      NaN    NaN   0;
                    2   200    0.05*200  200    NaN      NaN      NaN    NaN   0;
                    2   40     0.1*40    40     NaN      NaN      NaN    NaN   0;
                    0   r      0         r      NaN      NaN      NaN    NaN   0;
                    0   t      0         t      NaN      NaN      NaN    NaN   0];

probdata.marg(1,2:4)    = probdata.marg(1,2:4)*1000^2;
probdata.marg(2:3,2:4)  = probdata.marg(2:3,2:4)/1000;
                                   
% Correlation matrix 
probdata.correlation = eye(length(probdata.marg));

probdata.transf_type = 3;
probdata.Ro_method   = 1;
probdata.flag_sens   = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'ANALYSISOPT'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
analysisopt.echo_flag            = 0;

analysisopt.multi_proc           = 1;        % 1: block_size g-calls sent simultaneously
                                             %    - gfunbasic.m is used and a vectorized version of gfundata.expression is available.
                                             %      The number of g-calls sent simultaneously (block_size) depends on the memory
                                             %      available on the computer running FERUM.
                                             %    - gfunxxx.m user-specific g-function is used and able to handle block_size computations
                                             %      sent simultaneously, on a cluster of PCs or any other multiprocessor computer platform.
                                             % 0: g-calls sent sequentially
analysisopt.block_size           = 10^4;     % Number of g-calls to be sent simultaneously

% FORM analysis options
analysisopt.i_max                = 1000;     % Maximum number of iterations allowed in the search algorithm
analysisopt.e1                   = 1e-12;    % Tolerance on how close design point is to limit-state surface
analysisopt.e2                   = 1e-12;    % Tolerance on how accurately the gradient points towards the origin
analysisopt.step_code            = 0;        % 0: step size by Armijo rule, otherwise: given value is the step size
analysisopt.Recorded_u           = 1;        % 0: u-vector not recorded at all iterations, 1: u-vector recorded at all iterations
analysisopt.Recorded_x           = 1;        % 0: x-vector not recorded at all iterations, 1: x-vector recorded at all iterations

% FORM, SORM analysis options
analysisopt.grad_flag            = 'ffd';    % 'ddm': direct differentiation, 'ffd': forward finite difference
analysisopt.ffdpara              = 1000;     % Parameter for computation of FFD estimates of gradients - Perturbation = stdv/analysisopt.ffdpara;
                                             % Recommended values: 1000 for basic limit-state functions, 50 for FE-based limit-state functions
analysisopt.ffdpara_thetag       = 1000;     % Parameter for computation of FFD estimates of dbeta_dthetag
                                             % perturbation = thetag/analysisopt.ffdpara_thetag if thetag ~= 0 or 1/analysisopt.ffdpara_thetag if thetag == 0;
                                             % Recommended values: 1000 for basic limit-state functions, 100 for FE-based limit-state functions

% Simulation analysis (MC,IS,DS,SS) and distribution analysis options
analysisopt.num_sim              = 100000;   % Number of samples (MC,IS), number of samples per subset step (SS) or number of directions (DS)
analysisopt.rand_generator       = 1;        % 0: default rand matlab function, 1: Mersenne Twister (to be preferred)

% Simulation analysis (MC, IS) and distribution analysis options
analysisopt.sim_point            = 'dspt';  % 'dspt': design point, 'origin': origin in standard normal space (simulation analysis)
analysisopt.stdv_sim             = 1;        % Standard deviation of sampling distribution in simulation analysis
                                             
% Simulation analysis (MC, IS)
analysisopt.target_cov           = 0.01;   % Target coefficient of variation for failure probability
analysisopt.lowRAM               = 0;        % 1: memory savings allowed, 0: no memory savings allowed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'GFUNDATA' (one structure per gfun)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Type of limit-state function evaluator:
% 'basic': the limit-state function is defined by means of an analytical expression or a Matlab m-function,
%          using gfundata(lsf).expression. The function gfun.m calls gfunbasic.m, which evaluates gfundata(lsf).expression.
% 'xxx':   the limit-state function evaluation requires a call to an external code.  The function gfun.m calls gfunxxx.m,
%          which evaluates gfundata(lsf).expression where gext variable is a result of the external code.
gfundata(1).evaluator  = 'basic';
gfundata(1).type       = 'expression';   % Do not change this field!
        
% Expression of the limit-state function:
% L = 10m %% WARNING!
% a = 2m
gfundata(1).expression = 'rusting_beam_resistance(fy, b0, h0, r, t)';

% gfundata(1).dgdq       = { '  R ';
%                            '  theta_R ';
%                            '-(mu.*s_gn + g)*2*10^2/8';
%                            '-theta_E.*s_gn*2*10^2/8';
%                            '-theta_E.*mu*2*10^2/8';
%                            '-theta_E.*2*10^2/8'; };
                       
gfundata(1).thetag     =  [];

gfundata(1).thetagname = { };


% Flag for computation of sensitivities w.r.t. thetag parameters of the limit-state function
% 1: all sensitivities assessed, 0: no sensitivities assessment
gfundata(1).flag_sens  = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'FEMODEL'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

femodel = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'RANDOMFIELD'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randomfield = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FORM ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function updates probdata and gfundata before any analysis (must be run only once)
[probdata, gfundata, analysisopt] = update_data(1, probdata, analysisopt, gfundata, []);

% This function completely determines and updates parameters, mean and standard deviation associated with the distribution of each random variable
probdata.marg = distribution_parameter(probdata.marg);

% FORM analysis %
[formresults, ~] = form(1, probdata, analysisopt, gfundata, femodel, randomfield);

beta        = formresults.beta;

% formresults
% formresults.Recorded_x(5,:)
% for testing using Pf, WARNING!
% beta        = normcdf(-beta);

end