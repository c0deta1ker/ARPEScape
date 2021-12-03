function [fitStr0, statStr] = arpes2model_solver_n_runs(arpesStr, modelStr, iparams, ibgrnd, solve_type, kx_lims, eb_lims, n_runs)
% [fitStr0, statStr] = arpes2model_solver_n_runs(arpesStr, modelStr, iparams, ibgrnd, solve_type, kx_lims, eb_lims, n_runs)
%   Function that runs arpes2model_solver() N amount of independent times until
%   convergence. The initial conditions are randomly sampled over the range
%   of the input uncertainty. This allows the uncertainty in the best fit 
%   parameters to be determined by looking at the variance of the 
%   converged fit parameters.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   arpesStr:       MATLAB data-structure that contains the initial ARPES data.
%   -   modelStr:       MATLAB data-structure that contains the initial model data.
%   -   iparams:      	3 cells {x0}{lb}{ub} with 1x5 array: the model fit parameters [K0,KFWHM,EFWHM,MFP,BOFF,INT]
%   -   ibgrnd:         3 cells {x0}{lb}{ub} with 1x3 vector of the background parameters: [MX,MY,C]
%   -   solve_type:     String of either "fminunc", "fmincon", "lsqnonlin", "simulannealbnd" that sets the optimisation routine to use
%   -   kx_lims:      	[1x2] vector of [minKx, maxKx], which defines the consistent kx fit window
%   -   eb_lims:        [1x2] vector of [minEb, maxEb], which defines the consistent eb fit window
%   -   n_runs:         scalar value for the total number of independent runs
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information
%   -   statStr:        MATLAB data-structure that contains all the statistical analysis

%% Default parameters
% - Default input parameters to use
if nargin < 8; n_runs = 10; end
if nargin < 7; eb_lims = 0.95*[min(modelStr.eb(:)), max(modelStr.eb(:))]; end
if nargin < 6; kx_lims = 0.95*[min(modelStr.kx(:)), max(modelStr.kx(:))]; end
if nargin < 5; solve_type = "fmincon"; end
if nargin < 4; ibgrnd{1} = [0, 0, 0]; ibgrnd{2} = ibgrnd{1} - 0.10; ibgrnd{3} = ibgrnd{1} + 0.10; end
if isempty(n_runs); n_runs = 10; end
if isempty(eb_lims); eb_lims = 0.95*[min(modelStr.eb(:)), max(modelStr.eb(:))]; end
if isempty(kx_lims); kx_lims = 0.95*[min(modelStr.kx(:)), max(modelStr.kx(:))]; end
if isempty(solve_type); solve_type = "fmincon"; end
if isempty(ibgrnd); ibgrnd{1} = [0, 0, 0]; ibgrnd{2} = ibgrnd{1}-0.10; ibgrnd{3} = ibgrnd{1}+0.10; end
% - Total number of parameters
nParams = length(iparams{1}) + length(ibgrnd{1});

%% - 1 - Running the fitting algorithm N independent time for statistical analysis
statStr         = struct();
statStr.n_runs	= 1:n_runs;
fit     = cell(n_runs,1);
params	= zeros(n_runs, nParams);
CHISQ 	= zeros(n_runs,1);
% -- Running the fitting algorithm for N independent times for analysis
for i = statStr.n_runs
    fprintf("Run %i / %i \n", i, n_runs);
    % --- Add perturbation on the initial conditions
    new_init_conds          = iparams;
    new_init_conds{1}       = iparams{2} + (iparams{3}-iparams{2}) .* rand(size(iparams{2}));
    % --- Executing the fit
    if i == 1;  fitStr{1} = arpes2model_solver(arpesStr, modelStr, iparams, ibgrnd, solve_type, kx_lims, eb_lims); 
    else;       fitStr{i} = arpes2model_solver(arpesStr, modelStr, new_init_conds, ibgrnd, solve_type, kx_lims, eb_lims);
    end
    % --- Storing the best fits for statistical analysis
    fit{i}       	= fitStr{i};
    % --- Storing the best fit parameters for statistical analysis
    params(i,:)     = fit{i}.params';
    CHISQ(i)        = fit{i}.CHISQ;
end
% -- Assigning variables
statStr.fit     = fit;
statStr.params  = params;
statStr.CHISQ 	= CHISQ;
% -- Storing the parameter labels
statStr.params_label = {'K0', 'KFWHM', 'EFWHM', 'MFP', 'BOFF', 'INT', 'MX', 'MY', 'C'};

%% - 2 - Storing the fit structure that has the smallest value of chi-squared
[~, nCHISQ] = min(statStr.CHISQ(:));
fitStr0     = statStr.fit{nCHISQ};

%% - 3 - Extracting the best estimate for each one of the parameters and their standard deviation
% -- CHISQ
statStr.CHISQmu 	= mean(statStr.CHISQ);
statStr.CHISQerr	= 3*std(statStr.CHISQ);
% -- K0
statStr.K0          = statStr.params(:,1)';
statStr.K0mu        = mean(statStr.K0);
statStr.K0err       = 3*std(statStr.K0);
% -- kx FWHM
statStr.KFWHM       = statStr.params(:,2)';
statStr.KFWHMmu     = mean(statStr.KFWHM);
statStr.KFWHMerr    = 3*std(statStr.KFWHM);
% -- eb FWHM
statStr.EFWHM       = statStr.params(:,3)';
statStr.EFWHMmu  	= mean(statStr.EFWHM);
statStr.EFWHMerr	= 3*std(statStr.EFWHM);
% -- MFP
statStr.MFP         = statStr.params(:,4)';
statStr.MFPmu       = mean(statStr.MFP);
statStr.MFPerr      = 3*std(statStr.MFP);
% -- BOFF
statStr.BOFF        = statStr.params(:,5)';
statStr.BOFFmu      = mean(statStr.BOFF);
statStr.BOFFerr     = 3*std(statStr.BOFF);
% -- INT
statStr.INT         = statStr.params(:,6)';
statStr.INTmu       = mean(statStr.INT);
statStr.INTerr      = 3*std(statStr.INT);
% -- MX
statStr.MX          = statStr.params(:,end-2)';
statStr.MXmu        = mean(statStr.MX);
statStr.MXerr       = 3*std(statStr.MX);
% -- MY
statStr.MY          = statStr.params(:,end-1)';
statStr.MYmu        = mean(statStr.MY);
statStr.MYerr       = 3*std(statStr.MY);
% -- C
statStr.C           = statStr.params(:,end)';
statStr.Cmu         = mean(statStr.C);
statStr.Cerr        = 3*std(statStr.C);

end