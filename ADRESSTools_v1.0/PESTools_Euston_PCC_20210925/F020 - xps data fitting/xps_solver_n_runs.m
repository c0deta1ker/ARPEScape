function [fitStr0, statStr] = xps_solver_n_runs(xpsStr, cTYPE, iPESCurves, bTYPE, iPESBgrnd, solve_type, n_runs)
% [fitStr0, statStr] = xps_solver_n_runs(xpsStr, cTYPE, iPESCurves, bTYPE, iPESBgrnd, solve_type, n_runs)
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
    iPESBgrnd{1} = [min(xpsStr.xdat(:))+0.2, max(xpsStr.xdat(:))-0.2, 1, 0.5, 0, 0];
    iPESBgrnd{2} = abs(0.0.*iPESBgrnd{1}); iPESBgrnd{2}(:,6) = iPESBgrnd{1}(:,6) - 0.05;
    iPESBgrnd{3} = abs(0.0.*iPESBgrnd{1}); iPESBgrnd{2}(:,6) = iPESBgrnd{1}(:,6) + 0.05;
end
if nargin < 4; bTYPE = "Poly"; end
if isempty(n_runs); n_runs = 10; end
if isempty(solve_type); solve_type = "fmincon"; end
if isempty(iPESBgrnd)
    iPESBgrnd{1} = [min(xpsStr.xdat(:))+0.2, max(xpsStr.xdat(:))-0.2, 1, 0.5, 0, 0];
    iPESBgrnd{2} = abs(0.0.*iPESBgrnd{1}); iPESBgrnd{2}(:,6) = iPESBgrnd{1}(:,6) - 0.05;
    iPESBgrnd{3} = abs(0.0.*iPESBgrnd{1}); iPESBgrnd{2}(:,6) = iPESBgrnd{1}(:,6) + 0.05;
end
if isempty(bTYPE); bTYPE = "Poly"; end
% - Total number of states to be fitted
n_states = length(cTYPE);

%% - 1 - Running the fitting algorithm N independent time for statistical analysis
statStr         = struct();
statStr.n_runs	= 1:n_runs;
fit     = cell(n_runs,1);
params	= cell(10,1);
CHISQ 	= zeros(n_runs,1);
% - 1 - Running the fitting algorithm for N independent times for analysis
for i = statStr.n_runs
    fprintf("Run %i / %i", i, n_runs);
    % --- Add perturbation on the initial conditions
    new_init_conds      = iPESCurves;
    new_init_conds{1}   = iPESCurves{2} + (iPESCurves{3}-iPESCurves{2}) .* rand(size(iPESCurves{2}));
    % --- Executing the fit
    if i == 1;  fitStr{1} = xps_solver(xpsStr, cTYPE, iPESCurves, bTYPE, iPESBgrnd, solve_type);
    else;       fitStr{i} = xps_solver(xpsStr, cTYPE, new_init_conds, bTYPE, iPESBgrnd, solve_type);
    end
    % --- Storing the best fits for statistical analysis
    fit{i}          = fitStr{i};
    % --- Storing the best fit parameters for statistical analysis
    params{1}(i,:)	= fit{i}.BE;
    params{2}(i,:)	= fit{i}.INT;
    params{3}(i,:)	= fit{i}.FWHM;
    params{4}(i,:)	= fit{i}.MR;
    params{5}(i,:) 	= fit{i}.LSE;
    params{6}(i,:)	= fit{i}.LSI;
    params{7}(i,:) 	= fit{i}.LSW;
    params{8}(i,:)	= fit{i}.ASY;
    params{9}(i,:)	= fit{i}.AREA0;
    params{10}(i,:)	= fit{i}.BGR;
    params{11}(i,:)	= fit{i}.CHISQ;
    CHISQ(i)        = fit{i}.CHISQ;
end
% -- Assigning variables
statStr.fit     = fit;
statStr.params  = params;
statStr.CHISQ 	= CHISQ;
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