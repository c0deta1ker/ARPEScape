function [fitStr0, statStr] = arpes2model_solver_p_prctiles(arpesStr, modelStr, iparams, ibgrnd, solve_type, kx_lims, eb_lims, n_runs, prc, p_runs)
% [fitStr0, statStr] = arpes2model_solver_p_prctiles(modelStr, arpesStr, iparams, ibgrnd, solve_type, n_runs, prc, p_runs)
%   Function that runs arpes2model_solver() N amount of independent times until
%   convergence. The initial conditions are randomly sampled over the range
%   of the input uncertainty. This allows the uncertainty in the best fit 
%   parameters to be determined by looking at the variance of the 
%   converged fit parameters. Once it runs N independent times, it then
%   isolates the percentile of optimal fits based on chi-squared between 
%   (0 < prc). The initial conditions and bounds are re-calculated as the
%   mean and 3 sigma variation of the isolated fit parameters. It then
%   re-runs the arpes2model_solver() N amount of independent times again,
%   and repeats this process for 'p_runs' amount of times.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   arpesStr:       MATLAB data-structure that contains the initial ARPES data.
%   -   modelStr:       MATLAB data-structure that contains the initial model data.
%   -   iparams:      	3 cells {x0}{lb}{ub} with 1x5 array: the modelStr fit parameters [K0,KFWHM,EFWHM,MFP,BOFF,INT]
%   -   ibgrnd:         3 cells {x0}{lb}{ub} with 1x3 vector of the background parameters: [MX,MY,C]
%   -   solve_type:     String of either "fminunc", "fmincon", "lsqnonlin", "simulannealbnd" that sets the optimisation routine to use
%   -   kx_lims:      	[1x2] vector of [minKx, maxKx], which defines the consistent kx fit window
%   -   eb_lims:        [1x2] vector of [minEb, maxEb], which defines the consistent eb fit window
%   -   n_runs:         scalar value for the total number of independent runs
%   -   prc:            scalar value between [0, 100] that thresholds all the fit parameters within a percentile, 0 < prc, of the chi-squared values. This allows you to isolate the best solutions for diagnostic testing.
%   -   p_runs:         scalar value of the total number of percentile squeezes to perform.
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information
%   -   statStr:        MATLAB data-structure that contains all the statistical analysis

%% Default parameters
% - Default input parameters to use
if nargin < 10; p_runs = 2; end
if nargin < 9; prc = 10; end
if nargin < 8; n_runs = 10; end
if nargin < 7; eb_lims = 0.95*[min(modelStr.eb(:)), max(modelStr.eb(:))]; end
if nargin < 6; kx_lims = 0.95*[min(modelStr.kx(:)), max(modelStr.kx(:))]; end
if nargin < 5; solve_type = "fmincon"; end
if nargin < 4; ibgrnd{1} = [0, 0, 0]; ibgrnd{2} = ibgrnd{1} - 0.10; ibgrnd{3} = ibgrnd{1} + 0.10; end
if isempty(p_runs); p_runs = 10; end
if isempty(prc); prc = 2; end
if isempty(n_runs); n_runs = 10; end
if isempty(eb_lims); eb_lims = 0.95*[min(modelStr.eb(:)), max(modelStr.eb(:))]; end
if isempty(kx_lims); kx_lims = 0.95*[min(modelStr.kx(:)), max(modelStr.kx(:))]; end
if isempty(solve_type); solve_type = "fmincon"; end
if isempty(ibgrnd); ibgrnd{1} = [0, 0, 0]; ibgrnd{2} = ibgrnd{1}-0.10; ibgrnd{3} = ibgrnd{1}+0.10; end

%% - 1 - INITIALISING VARIABLES
statStr = {};

%% - 2 - RUNNING RUNNING THROUGH THE PERCENTILE-CONVERGING FITTING METHOD
for p = 1:p_runs
    % -- Using the arpes2model solver for n runs
    [fitStr0, statStr{p}] = arpes2model_solver_n_runs(arpesStr, modelStr, iparams, ibgrnd, solve_type, kx_lims, eb_lims, n_runs);
    % -- Extracting the percentile value/index for the best solutions
    chisq_prctile   = prctile(statStr{p}.CHISQ, prc);
    prctile_indx    = statStr{p}.CHISQ<chisq_prctile;
    if isempty(prctile_indx) || sum(prctile_indx) == 0 || range(chisq_prctile) == 0 || isempty(range(chisq_prctile))
        [~, prctile_indx]   = min(statStr{p}.CHISQ(:));
    end
    fitStr0.iter        = p;
    fitStr0.prc_exit    = 0;
    % -- If only a single value of chi-squared remains, this is the best solution
    if length(prctile_indx) == 1
        fitStr0.prc_exit = 1;
        break
    % -- If there are multiple values, redefine the parameters and run again 
    else
        % -- Recalculating the new best fit parameters and bounds
        params      = statStr{p}.params(prctile_indx,:);
        mu_params   = mean(params, 1);
        std_params  = 3.*std(params, 1);
        % -- Defining the new best fit parameters
        iparams     = {};
        iparams{1}  = mu_params(1:end-3);
        iparams{2}  = mu_params(1:end-3) - std_params(1:end-3);
        iparams{3}  = mu_params(1:end-3) + std_params(1:end-3);
        % -- Defining the new best fit background parameters
        ibgrnd      = {};
        ibgrnd{1}   = mu_params(end-2:end);
        ibgrnd{2}   = mu_params(end-2:end) - std_params(end-2:end);
        ibgrnd{3}   = mu_params(end-2:end) + std_params(end-2:end);
        % -- Verifying that the bound limits are not identical
        for i = 1:size(iparams{1}, 2)
            if iparams{1}(i) == iparams{2}(i); iparams{2}(i) = iparams{1}(i) - 0.001; end
            if iparams{1}(i) == iparams{3}(i); iparams{3}(i) = iparams{1}(i) + 0.001; end
        end
        for i = 1:size(ibgrnd{1}, 2)
            if ibgrnd{1}(i) == ibgrnd{2}(i); ibgrnd{2}(i) = ibgrnd{1}(i) - 0.001; end
            if ibgrnd{1}(i) == ibgrnd{3}(i); ibgrnd{3}(i) = ibgrnd{1}(i) + 0.001; end
        end
%         % -- Verifying that the bound limits are not too small
%         for i = 1:size(iparams{1}, 2)
%             if abs(iparams{1}(i) - iparams{2}(i)) < 0.001; iparams{2}(i) = iparams{1}(i) - 0.001; end
%             if abs(iparams{1}(i) - iparams{3}(i)) < 0.001; iparams{3}(i) = iparams{1}(i) - 0.001; end
%         end
%         for i = 1:size(ibgrnd{1}, 2)
%             if abs(ibgrnd{1}(i) - ibgrnd{2}(i)) < 0.001; ibgrnd{2}(i) = ibgrnd{1}(i) - 0.001; end
%             if abs(ibgrnd{1}(i) - ibgrnd{3}(i)) < 0.001; ibgrnd{3}(i) = ibgrnd{1}(i) - 0.001; end
%         end
    end
end
end