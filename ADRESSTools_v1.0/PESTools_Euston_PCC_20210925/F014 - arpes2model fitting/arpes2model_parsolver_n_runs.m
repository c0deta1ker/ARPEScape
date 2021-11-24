function [fitStr0, statStr] = arpes2model_parsolver_n_runs(MODEL, ARPES, iparams, ibgrnd, solve_type, n_runs)
% [fitStr0, statStr] = arpes2model_solver_n_runs(MODEL, ARPES, iparams, ibgrnd, solve_type, n_runs)
%   Function that runs arpes2model_solver() N amount of independent times until
%   convergence. The initial conditions are randomly sampled over the range
%   of the input uncertainty. This allows the uncertainty in the best fit 
%   parameters to be determined by looking at the variance of the 
%   converged fit parameters.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   MODEL:          MATLAB data-structure that contains the model data.
%   -   ARPES:         	MATLAB data-structure that contains the ARPES data.
%   -   iparams:      	3 cells {x0}{lb}{ub} with 1x5 array: the model fit parameters [kxFWHM,ebFWHM,MFP,BOFF,INT]
%   -   ibgrnd:         3 cells {x0}{lb}{ub} with 1x3 vector of the background parameters: [MX,MY,C]
%   -   solve_type:     String of either "lsqnonlin", "simulannealbnd"
%   -   n_runs:         scalar value for the total number of independent runs
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information
%   -   statStr:        MATLAB data-structure that contains all the statistical analysis

%% Default parameters
% - Default input parameters to use
if nargin < 6; n_runs = 10; end
if nargin < 5; solve_type = "fmincon"; end
if nargin < 4; ibgrnd = [0, 0, 0]; end
if isempty(n_runs); n_runs = 10; end
if isempty(solve_type); solve_type = "fmincon"; end
if isempty(ibgrnd); ibgrnd = [0, 0, 0]; end
% - Total number of parameters
nParams = length(iparams{1}) + length(ibgrnd{1});

%% - 1 - Running the fitting algorithm N independent time for statistical analysis
% -- Initialising parallel processing
nCPU_PCC_Work = 6;
nCPU_PSI_HPC_Merlin6 = 44;
myCluster = parcluster('local'); nCPU = myCluster.NumWorkers;
if isempty(gcp('nocreate')); poolobj = parpool('local', nCPU, 'IdleTimeout', Inf); end
% -- Initialising variables
statStr         = struct();
statStr.n_runs	= 1:n_runs;
fit     = cell(n_runs,1);
params	= zeros(n_runs, nParams);
CHISQ 	= zeros(n_runs,1);
% -- Running the fitting algorithm N independent times for analysis
parfor i = statStr.n_runs
    if n_runs < 101; fprintf("Run %i / %i \n", i, n_runs); end
    % --- Uniform sampling of the initial conditions
    new_init_conds          = iparams;
    new_init_conds{1}       = iparams{2} + (iparams{3}-iparams{2}) .* rand(size(iparams{2}));
    % --- Executing the fit
    fitStr = arpes2model_solver(MODEL, ARPES, new_init_conds, ibgrnd, solve_type);
    % --- Storing the best fit parameters for statistical analysis
    fit{i}       	= fitStr;
    % --- Storing the best fit parameters for statistical analysis
    params(i,:)     = fit{i}.params';
    CHISQ(i)        = fit{i}.CHISQ;
end
% -- Shutting down parallel processing
% delete(poolobj);
% -- Assigning variables
statStr.fit     = fit;
statStr.params  = params;
statStr.CHISQ 	= CHISQ;
% -- Storing the parameter labels
statStr.params_label = {'kxFWHM', 'ebFWHM', 'MFP', 'BOFF', 'INT', 'MX', 'MY', 'C'};

%% - 2 - Storing the fit structure that has the smallest value of chi-squared
[~, nCHISQ] = min(statStr.CHISQ(:));
fitStr0     = statStr.fit{nCHISQ};

%% - 3 - Extracting the best estimate for each one of the parameters and their standard deviation
% -- CHISQ
statStr.CHISQmu 	= mean(statStr.CHISQ);
statStr.CHISQerr	= 3*std(statStr.CHISQ);
% -- kx FWHM
statStr.kxFWHM      = statStr.params(:,1)';
statStr.kxFWHMmu    = mean(statStr.kxFWHM);
statStr.kxFWHMerr   = 3*std(statStr.kxFWHM);
% -- eb FWHM
statStr.ebFWHM      = statStr.params(:,2)';
statStr.ebFWHMmu  	= mean(statStr.ebFWHM);
statStr.ebFWHMerr	= 3*std(statStr.ebFWHM);
% -- MFP
statStr.MFP         = statStr.params(:,3)';
statStr.MFPmu       = mean(statStr.MFP);
statStr.MFPerr      = 3*std(statStr.MFP);
% -- BOFF
statStr.BOFF        = statStr.params(:,4)';
statStr.BOFFmu      = mean(statStr.BOFF);
statStr.BOFFerr     = 3*std(statStr.BOFF);
% -- INT
statStr.INT         = statStr.params(:,5)';
statStr.INTmu       = mean(statStr.INT);
statStr.INTerr      = 3*std(statStr.INT);
% -- MX
statStr.MX          = statStr.params(:,end-3)';
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