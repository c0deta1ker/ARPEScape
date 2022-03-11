function [fitStr0, statStr] = pes2ncurve_parsolver_n_runs(pesStr, cTYPE, iparams, bTYPE, ibgrnd, solve_type, n_runs)
% [fitStr0, statStr] = pes2ncurve_parsolver_n_runs(xpsStr, cTYPE, iPESCurves, bTYPE, iPESBgrnd, solve_type, n_runs)
%   Function that runs 'xps_solver()' N amount of independent times until
%   convergence. The initial conditions are randomly sampled over the range
%   of the input uncertainty. This allows the uncertainty in the best fit 
%   parameters to be determined by looking at the variance of the 
%   converged fit parameters. It is advised to run this algorithm for AT
%   LEAST N = 100, so good statistics can be obtained. Ideally, N = 1000
%   is optimal, but this can take a few hours, especially when using
%   the 'simulannealbnd()' global optimisation method.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xpsStr:         MATLAB data-structure that contains the XPS data.
%   -   cTYPE:          1xN vector of the type of curve to use for the nth state.
%   -   iPESCurves:   	3 cells {x0}{lb}{ub} with Nx8 arrays: the n'th peak parameters [BE,INT,FWHM,MR,LSE,LSI,LSW,ASY]
%   -   bTYPE:          string of the type of background to use for fitting. Default: "Poly" ("none", "Shir", "LinShir")
%   -   iPESBgrnd:      3 cells {x0}{lb}{ub} with 1x5 vectors: the background parameters: [LHS,RHS,ORD,LAM,DEL,BGR]
%   -   solve_type:     string of the optimisation algorithm to use; "lsqcurvefit" (very fast, good), "lsqnonlin" (fast, not good), "simulannealbnd" (slow, best). 
%   -   n_runs:         scalar value for the total number of independent runs
%
%   OUT:
%   -   fitStr0:        MATLAB data-structure that contains all the fit parameters / variables / information for the fit with minimum CHISQ value
%   -   statStr:        MATLAB data-structure that contains all the statistical analysis

%% Default parameters
% - Default input parameters to use
if nargin < 7; n_runs = 10; end
if nargin < 6; solve_type = "fmincon"; end
if nargin < 5
    ibgrnd{1} = [min(pesStr.xdat(:))+0.2, max(pesStr.xdat(:))-0.2, 1, 0.5, 0, 0];
    ibgrnd{2} = abs(0.0.*ibgrnd{1}); ibgrnd{2}(:,6) = ibgrnd{1}(:,6) - 0.05;
    ibgrnd{3} = abs(0.0.*ibgrnd{1}); ibgrnd{2}(:,6) = ibgrnd{1}(:,6) + 0.05;
end
if nargin < 4; bTYPE = "Poly"; end
if isempty(n_runs); n_runs = 10; end
if isempty(solve_type); solve_type = "fmincon"; end
if isempty(ibgrnd)
    ibgrnd{1} = [min(pesStr.xdat(:))+0.2, max(pesStr.xdat(:))-0.2, 1, 0.5, 0, 0];
    ibgrnd{2} = abs(0.0.*ibgrnd{1}); ibgrnd{2}(:,6) = ibgrnd{1}(:,6) - 0.05;
    ibgrnd{3} = abs(0.0.*ibgrnd{1}); ibgrnd{2}(:,6) = ibgrnd{1}(:,6) + 0.05;
end
if isempty(bTYPE); bTYPE = "Poly"; end
% - Total number of states to be fitted
n_states = length(cTYPE);

%% - 1 - Running the fitting algorithm N independent time for statistical analysis
% -- Initialising parallel processing
nCPU_PCC_Work = 6;
nCPU_PSI_HPC_Merlin6 = 44;
myCluster = parcluster('local'); nCPU = myCluster.NumWorkers;
if isempty(gcp('nocreate')); poolobj = parpool('local', nCPU, 'IdleTimeout', Inf); end
% -- Initialising data structure
statStr         = struct();
statStr.n_runs	= 1:n_runs;
% -- Initialising variables
fit     = cell(n_runs,1);
BE      = zeros(n_runs,n_states);
INT     = zeros(n_runs,n_states);
FWHM	= zeros(n_runs,n_states);
MR      = zeros(n_runs,n_states);
LSE 	= zeros(n_runs,n_states);
LSI     = zeros(n_runs,n_states);
LSW 	= zeros(n_runs,n_states);
ASY     = zeros(n_runs,n_states);
AREA0	= zeros(n_runs,n_states);
BGR     = zeros(n_runs,1);
CHISQ 	= zeros(n_runs,1);
% - 1 - Running the fitting algorithm for N independent times for analysis
parfor i = statStr.n_runs
    fprintf("Run %i / %i", i, n_runs);
    % --- Add perturbation on the initial conditions
    new_init_conds      = iparams;
    new_init_conds{1}   = iparams{2} + (iparams{3}-iparams{2}) .* rand(size(iparams{2}));
    % --- Executing the fit
    fitStr{i} = xps_solver(pesStr, cTYPE, new_init_conds, bTYPE, ibgrnd, solve_type);
    % --- Storing the best fits for statistical analysis
    fit{i}          = fitStr{i};
    % --- Storing the best fit parameters for statistical analysis
    BE(i,:)     = fit{i}.BE;
    INT(i,:)	= fit{i}.INT;
    FWHM(i,:)	= fit{i}.FWHM;
    MR(i,:)     = fit{i}.MR;
    LSE(i,:) 	= fit{i}.LSE;
    LSI(i,:)	= fit{i}.LSI;
    LSW(i,:) 	= fit{i}.LSW;
    ASY(i,:)	= fit{i}.ASY;
    AREA0(i,:)	= fit{i}.AREA0;
    BGR(i)      = fit{i}.BGR;
    CHISQ(i)  	= fit{i}.CHISQ;
end
% -- Shutting down parallel processing
% delete(poolobj);
% -- Assigning variables
statStr.fit     = fit;
% --- Storing the best fit parameters for statistical analysis
statStr.params      = {};
statStr.params{1}	= BE;
statStr.params{2}	= INT;
statStr.params{3}	= FWHM;
statStr.params{4}	= MR;
statStr.params{5} 	= LSE;
statStr.params{6}	= LSI;
statStr.params{7} 	= LSW;
statStr.params{8}	= ASY;
statStr.params{9}	= AREA0;
statStr.params{10}	= BGR;
statStr.params{11}	= CHISQ;
statStr.CHISQ       = CHISQ;
% -- Storing the parameter labels
statStr.params_label = {'BE', 'INT', 'FWHM', 'MR', 'LSE', 'LSI', 'LSW', 'ASY', 'AREA0', 'BGR', 'CHISQ'};

%% - 2 - Storing the fit structure that has the smallest value of chi-squared
[~, nCHISQ] = min(statStr.CHISQ(:));
fitStr0     = statStr.fit{nCHISQ};

%% - 3 - Extracting the best estimate for each one of the parameters and their standard deviation
% -- CHISQ
statStr.CHISQmu   	= mean(statStr.CHISQ);
statStr.CHISQerr  	= 3*std(statStr.CHISQ);
% -- Binding energy
statStr.BE          = statStr.params{1};
statStr.BEmu        = mean(statStr.BE, 1);
statStr.BEerr       = 3*std(statStr.BE, 1);
% -- Intensity
statStr.INT         = statStr.params{2};
statStr.INTmu       = mean(statStr.INT, 1);
statStr.INTerr      = 3*std(statStr.INT, 1);
% -- FWHM
statStr.FWHM        = statStr.params{3};
statStr.FWHMmu      = mean(statStr.FWHM, 1);
statStr.FWHMerr   	= 3*std(statStr.FWHM, 1);
% -- MR
statStr.MR          = statStr.params{4};
statStr.MRmu        = mean(statStr.MR, 1);
statStr.MRerr       = 3*std(statStr.MR, 1);
% -- LSE
statStr.LSE         = statStr.params{5};
statStr.LSEmu       = mean(statStr.LSE, 1);
statStr.LSEerr      = 3*std(statStr.LSE, 1);
% -- LSI
statStr.LSI         = statStr.params{6};
statStr.LSImu       = mean(statStr.LSI, 1);
statStr.LSIerr      = 3*std(statStr.LSI, 1);
% -- LSW
statStr.LSW         = statStr.params{7};
statStr.LSWmu       = mean(statStr.LSW, 1);
statStr.LSWerr      = 3*std(statStr.LSW, 1);
% -- ASY
statStr.ASY         = statStr.params{8};
statStr.ASYmu       = mean(statStr.ASY, 1);
statStr.ASYerr      = 3*std(statStr.ASY, 1);
% -- AREA0
statStr.AREA0       = statStr.params{9};
statStr.AREA0mu     = mean(statStr.AREA0, 1);
statStr.AREA0err    = 3*std(statStr.AREA0, 1);
% -- BGR
statStr.BGR         = statStr.params{10};
statStr.BGRmu       = mean(statStr.BGR, 1);
statStr.BGRerr      = 3*std(statStr.BGR, 1);

end