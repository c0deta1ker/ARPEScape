function [fitStr, Ef] = fdd2fit_solver(Angle, Energy, Data, AngleWin, EnergyWin, fddType, iparams, solve_type, plot_result)
% [fitStr, Ef] = fdd2fit_solver(Angle, Energy, Data, fddType, iparams, solve_type)
%   Function that runs an optimisation algorithm to solve for the best fit
%   to a Fermi-Dirac Distribution. The FDD function can be chosen for
%   optimal fitting to the data.
%   
%   REQ. FUNCTIONS: none
%   
%   IN:
%   -   Angle:              [1×nAngle] row vector of angles.
%   -   Energy:             [nEnergyx1] column vector of energy.
%   -   Data:               [nEnergy x nAngle] data array.
%   -   AngleWin:           1x2 vector of [minAngle, maxAngle] over which the fit is performed. Default: 95% of the full Angle window.
%   -   EnergyWin:          1x2 vector of [minEnergy, maxEnergy] over which the fit is performed. Default: 95% of the full Energy window.
%   -   fddType:            string of the type of FDD curve to use for fitting. Default: "FDDGsL" ("FDD", "FDDG", "FDDGpL", "FDDGsL")
%   -   iparams:            3 cells {x0}{lb}{ub} with 1×7 vectors: defines the input parameters [fdd_ef,fdd_T,fdd_fwhm,lin_grad,lin_offset,const_bgrnd,fdd_scale]
%   -   solve_type:         string of the optimisation algorithm to use; "lsqcurvefit" (very fast, good), "lsqnonlin" (fast, not good), "simulannealbnd" (slow, best). 
%   -   plot_result:        if True, will plot the result of the fitting.
%
%   OUT:
%   -   fitStr:             MATLAB data-structure that contains all the fit parameters / variables / information
%   -   Ef:                 scalar value of the position of the Fermi-level.

%% Default parameters
if nargin < 9; plot_result = 1; end
if isempty(plot_result); plot_result = 1; end
% - Solver parameters
if nargin < 8; solve_type = "lsqcurvefit"; end
if nargin < 6; fddType = "FDDGsL"; end
if isempty(solve_type); solve_type = "lsqcurvefit"; end
if isempty(fddType); fddType = "FDDGsL"; end
% -- Window parameters
if ~isempty(Angle)
    if nargin < 4; AngleWin = mean(Angle(:)) + abs(0.475*range(Angle(:))).*[-1, 1]; end
    if isempty(AngleWin); AngleWin = mean(Angle(:)) + abs(0.475*range(Angle(:))).*[-1, 1]; end
else; AngleWin = [];
end
if nargin < 5; EnergyWin = mean(Energy(:)) + abs(0.475*range(Energy(:))).*[-1, 1]; end
if isempty(EnergyWin); EnergyWin = mean(Energy(:)) + abs(0.475*range(Energy(:))).*[-1, 1]; end
% -- Fit variables
if nargin < 7 || isempty(iparams)
    iparams{1} = [mean(EnergyWin(:)), 12,  0.30,  0.00, 0,               mean(Data(:)), 1];
    iparams{2} = [min(EnergyWin(:)),  12,  0.00, -5.00, min(Energy(:)), -max(Data(:)), 0];
    iparams{3} = [max(EnergyWin(:)),  12,  2.00,  0.00, max(Energy(:)),  max(Data(:)), max(Data(:))];
end
% -- Consistency check and finding the total number of curves
if size(iparams{1}, 1) ~= size(iparams{2}, 1) || size(iparams{2}, 1) ~= size(iparams{3}, 1)
    error('The input parameter cell array is not a consistent size - check iparams input!');
end

%% - 1 - CONSTRAINING THE FIT PARAMETERS TO PHYSICAL VALUES
% -- Defining the initial conditions and limits of the optimisation method
x0 = [iparams{1}(:)];
lb = [iparams{2}(:)];
ub = [iparams{3}(:)];
% -- Defining the initial conditions and limits of the optimisation method
% --- lower bounds
if lb(1) < min(EnergyWin(:)); lb(1) = min(EnergyWin(:)); end    % if fdd_ef < 0, then make it 0
if lb(2) < 0; lb(2) = 0; end    % if fdd_T < 0, then make it 0
if lb(3) < 0; lb(3) = 0; end    % if fdd_fwhm < 0, then make it 0
if lb(7) < 0; lb(7) = 0; end    % if fdd_scale < 0, then make it 0
% --- upper bounds
if ub(1) > max(EnergyWin(:)); ub(1) = max(EnergyWin(:)); end    % if fdd_ef < 0, then make it 0
if ub(4) > 0; ub(4) = 0; end        % if lin_grad > 0, then make it 0
% --- parameters
if x0(1) < lb(1); x0(1) = lb(1); elseif x0(1) > ub(1); x0(1) = ub(1); end
if x0(2) < lb(2); x0(2) = lb(2); elseif x0(2) > ub(2); x0(2) = ub(2); end
if x0(3) < lb(3); x0(3) = lb(3); elseif x0(3) > ub(3); x0(3) = ub(3); end
if x0(4) < lb(4); x0(4) = lb(4); elseif x0(4) > ub(4); x0(4) = ub(4); end
if x0(5) < lb(5); x0(5) = lb(5); elseif x0(5) > ub(5); x0(5) = ub(5); end
if x0(6) < lb(6); x0(6) = lb(6); elseif x0(6) > ub(6); x0(6) = ub(6); end
if x0(7) < lb(7); x0(7) = lb(7); elseif x0(7) > ub(7); x0(7) = ub(7); end
% -- Determining the degrees of freedom for the fitting model
val	= lb + (ub-lb); DoF = 0;
for i = 1:length(val); if val(i) ~= lb(i); DoF = DoF + 1; end; end

%% - 2 - EXTRACTING THE FDD DATA
% -- Cross-correlating sum if necessary
if size(Data,3) > 1; [Data,~]      = SumScanXC(Energy, Data); end
% -- Angle-Integrated sum if necessary
if ~isempty(Angle);         [ydat, xdat] = IntAngle(Data, Angle, Energy, AngleWin);
else;                       xdat = Energy; ydat = Data;
end
% -- Ensuring data is in a column vector
if size(xdat, 2) > 1; xdat = xdat'; end
if size(ydat, 2) > 1; ydat = ydat'; end
% -- Normalising the data
ydat(isnan(ydat)) = 0;
ydat = ydat - min(ydat(:)); 
ydat = ydat ./ max(ydat(:));
% -- Extracting the data to be fitted
[roi_xdat, roi_ydat] = data_crop1D(xdat, ydat, EnergyWin);
% -- Defining a structure that stores all relevant model and data variables
FDDObj              = struct();
FDDObj.xdat         = xdat;
FDDObj.ydat   	    = ydat;
FDDObj.roi_xdat     = roi_xdat;
FDDObj.roi_ydat   	= roi_ydat;
% -- Appending the input arguments to the global variable
FDDObj.fit_args.solve_type  = solve_type;
FDDObj.fit_args.fddType    	= fddType;
FDDObj.fit_args.iparams     = iparams;
FDDObj.fit_args.AngleWin    = AngleWin;
FDDObj.fit_args.EnergyWin   = EnergyWin;

%% - 3 - RUNNING THE SOLVER
% -- Local solver properties
MaxFunEvals = 1e5;
MaxIter     = 5e3;
FuncTol     = 1e-7; % 1e-7;
StepTol     = 1e-7; % 1e-7;
OptTol      = 1e-7; % 1e-7;
ConTol      = 1e-7; % 1e-7;
FinDiffRelStep = 1e-2; % 1e-5;
% -- Simulated annealing properties
TempFcn     = 'temperaturefast';
InitTemp    = 1e2;
ReannealInt = 1e2;
%% - 3.1 - LOCAL SOLVER: BOUNDED LEAST SQUARES FITTING METHOD FOR NON-LINEAR CURVE FITTING
if solve_type == "lsqcurvefit"
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions(@lsqcurvefit,...
        'MaxFunctionEvaluations', MaxFunEvals,...
        'FunctionTolerance', FuncTol,...
        'StepTolerance', StepTol,...
        'OptimalityTolerance', OptTol,...
        'MaxIterations', MaxIter,...
        'FinDiffRelStep', FinDiffRelStep,...
        'FinDiffType', 'central');
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,rnorm,resid,exitflag,output,lambda,jacobian]  = lsqcurvefit(@(x,roi_xdat) full_fdd_function(x,roi_xdat,FDDObj), x0, roi_xdat, roi_ydat, lb, ub, options);  

%% - 3.2 - LOCAL SOLVER: UNBOUNDED LEAST SQUARES NONLINEAR REGRESSION FIT
elseif solve_type == "nlinfit"
    % -- Defining the optimisation options for the simulated annealing method
    options = statset(...
        'MaxFunEvals', MaxFunEvals,...
        'MaxIter', MaxIter,...
        'TolFun', FuncTol);
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(roi_xdat, roi_ydat, @(x,roi_xdat) full_fdd_function(x,roi_xdat,FDDObj), x0, options);
%% - 3.3 - LOCAL SOLVER: FIND MINIMUM OF UNBOUNDED MULTIVARIATE FUNCTION
elseif solve_type == "fminunc"
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions('fminunc',...
        'MaxFunctionEvaluations', MaxFunEvals,...
        'FunctionTolerance', FuncTol,...
        'StepTolerance', StepTol,...
        'OptimalityTolerance', OptTol,...
        'MaxIterations', MaxIter,...
        'FinDiffRelStep', FinDiffRelStep,...
        'FinDiffType', 'central');
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,fval,exitflag,output,grad,hessian] = fminunc(@(x) minimize_function(x,FDDObj), x0, options); 

%% - 3.4 - LOCAL SOLVER: FIND MINIMUM OF BOUNDED MULTIVARIATE FUNCTION
elseif solve_type == "fmincon"
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions('fmincon',...
        'MaxFunctionEvaluations', MaxFunEvals,...
        'FunctionTolerance', FuncTol,...
        'StepTolerance', StepTol,...
        'OptimalityTolerance', OptTol,...
        'MaxIterations', MaxIter,...
        'ConstraintTolerance', ConTol,...
        'FinDiffRelStep', FinDiffRelStep,...
        'FinDiffType', 'central');
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) minimize_function(x,FDDObj), x0, [], [], [], [], lb, ub, [], options);  

    %% - 3.5 - LOCAL SOLVER: BOUNDED LEAST SQUARES FITTING METHOD FOR NON-LINEAR LEAST SQUARES PROBLEMS
elseif solve_type == "lsqnonlin"
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions(@lsqnonlin,...
        'Algorithm', 'levenberg-marquardt',...
        'MaxFunctionEvaluations', MaxFunEvals,...
        'FunctionTolerance', FuncTol,...
        'StepTolerance', StepTol,...
        'OptimalityTolerance', OptTol,...
        'MaxIterations', MaxIter);
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,rnorm,resid,exitflag,output,lambda,jacobian] = lsqnonlin(@(x) minimize_function(x,FDDObj), x0, lb, ub, options);  
%% - 3.8 - GLOBAL SOLVER: BOUNDED SIMULATED ANNEALING TO MINIMISE RESIDUALS
elseif solve_type == "simulannealbnd"
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    hybridopts = optimoptions('fmincon',...
        'MaxFunctionEvaluations', MaxFunEvals,...
        'OptimalityTolerance', OptTol,...
        'StepTolerance', StepTol,...
        'ConstraintTolerance', ConTol);
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions(@simulannealbnd,...
        'TemperatureFcn', TempFcn,...
        'InitialTemperature', InitTemp,...
        'ReannealInterval', ReannealInt,...
        'MaxFunctionEvaluations', MaxFunEvals,...
        'MaxIterations', MaxIter,...
        'TolFun', FuncTol,...
        'FunctionTolerance', FuncTol,...
        'HybridFcn' , {'fmincon' , hybridopts});
    % to view annealing, add:
    % 'PlotFcn',{@saplotbestf,@saplottemperature,@saplotf,@saplotstopping}, 'MaxIterations', 50,
    % -- Run the simulated annealing method
    % --- (1) Initially optimise all parameters
    [params,fval,exitflag,output] = simulannealbnd(@(x) minimize_function(x,FDDObj), x0, lb, ub, options);
end

%% - 4 - STORING THE FITS AS A MATLAB STRUCTURE
fitStr              = struct();
%% 4.1 - Storing the initial arguments
fitStr.solve_type	= solve_type;
fitStr.fddType      = fddType;
fitStr.iparams      = iparams;
fitStr.AngleWin     = AngleWin;
fitStr.EnergyWin    = EnergyWin;
%% 4.2 - Storing the XPS data, background, model and minimisation variables
% -- Extracting the original data 
fitStr.xdat     = xdat;
fitStr.ydat 	= ydat;
% -- Extracting the X-DOMAIN and DATA
[fitStr.X, fitStr.D] = fit_data(params, FDDObj);
% -- Extracting the MODEL
fitStr.M        = fit_model(params, FDDObj);
% -- Extracting the RESIDUALS (R = D - (M + B))
fitStr.R        = fitStr.D - fitStr.M;
% -- Extracting the CHI-SQUARED value to be minimised
fitStr.CHISQ    = sum(fitStr.R.^2 ./ abs(fitStr.M));
fitStr.MINFUN   = minimize_function(params, FDDObj);
% -- Append the degrees of freedom
fitStr.DoF   	= DoF;
%% 4.3 - Storing the final fit variables of each PES curve component
fitStr.XX           = linspace(min(fitStr.X(:))-1, max(fitStr.X(:))+1, 1e3)';
if fitStr.fddType == "FDD"; fitStr.YY = params(7).*FDD(fitStr.XX, params(1), params(2));
elseif fitStr.fddType == "FDDG"; fitStr.YY = params(7).*FDDG(fitStr.XX, params(1), params(2), params(3), params(4));
elseif fitStr.fddType == "FDDGpL"; fitStr.YY = params(7).*FDDGpL(fitStr.XX, params(1), params(2), params(3), params(4), params(5), params(6));
elseif fitStr.fddType == "FDDGsL"; fitStr.YY = params(7).*FDDGsL(fitStr.XX, params(1), params(2), params(3), params(4), params(5), params(6));
end
% --- Storing the best fit parameters
fitStr.params   = params;
fitStr.EF       = params(1); Ef = fitStr.EF;
fitStr.TEMP     = params(2);
fitStr.FWHM     = params(3);
fitStr.GRAD     = params(4);
fitStr.OFFSET   = params(5);
fitStr.BGRND    = params(6);
fitStr.SCALE    = params(7);

%% - 5 - PLOTTING THE FIT IF REQUIRED 
if plot_result == 1
    fig = figure(); 
    fig.Position(3) = 1.0*400;
    fig.Position(4) = 1.0*400;
    hold on;
    % -- Plot the figure
    plot(fitStr.xdat, fitStr.ydat, 'o', 'Color', [0.5 0.5 0.5], 'markerfacecolor', [0.5 0.5 0.5]);
    plot(fitStr.X, fitStr.D,'o', 'Color', 'r', 'markerfacecolor', 'r');
    plot(fitStr.XX, fitStr.YY, 'k-', 'linewidth', 3.0, 'linestyle', '-');
    line([1, 1].*fitStr.EF, [-1e5, 1e5], 'color', 'b', 'linewidth', 1.5, 'linestyle', '-');
    a = line([0 0], [-1e5, 1e5], 'Color', [0 0 0], 'LineWidth', 1, 'Linestyle', '--');
    b = line([-1e5, 1e5], [0 0], 'Color', [0 0 0], 'LineWidth', 1, 'Linestyle', '--');
    a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % -- Formatting the figure
    gca_props(0);
    title('fdd2fit_solver()', 'interpreter', 'none'); 
    ylabel('$$ \bf Intensity\ [arb.] $$', 'Interpreter', 'latex');
    xlabel('$$ \bf E_B\ [eV] $$', 'Interpreter', 'latex');
    legend({'', 'Data', 'Fit', 'Ef'}, 'location', 'northeast', 'interpreter', 'none');
    axis([...
        min(fitStr.X(:)), ...
        max(fitStr.X(:)), ...
        min(fitStr.D(:)), ...
        max(fitStr.D(:))]);
    % -- Add annotation for the quality of fit
    text(0.04, 0.35, "$$ \chi^2 = $$ " + string(fitStr.CHISQ), 'interpreter', 'latex', 'fontsize', 12, 'color', 'k', 'Units','normalized');
    text(0.04, 0.25, "$$ E_F = $$ " + string(round(fitStr.EF,3)) + " eV", 'interpreter', 'latex', 'fontsize', 12, 'color', 'k', 'Units','normalized');
   % -- Printing data structure
   Ef
end
end

%% DEFINING THE RESIDUAL MINIMISATION FUNCTION TO BE MINIMISED
function MINFUN = minimize_function(x, FDDObj)
    % - 1 - Extracting the DATA
    [~, D] = fit_data(x, FDDObj);
    % - 2 - Extracting the MODEL
    M = fit_model(x, FDDObj);
    % - 3 - Extracting the RESIDUALS (residuals = data - (model + background))
    R   = D - M;
    % - 4 - Extracting the CHI-SQUARED value to be minimised
    MINFUN = sum(R.^2 ./ abs(M));
end

%% DEFINING THE FUNCTION THAT EXTRACTS THE FDD DATA TO BE FITTED
function [X, D] = fit_data(x, FDDObj)
    X = FDDObj.roi_xdat;
    D = FDDObj.roi_ydat;
    D(isnan(D)) = 0; if size(D, 2) > 1; D = D'; end
end

%% DEFINING THE FUNCTION THAT DETERMINES THE FDD CURVE FIT
function M = fit_model(x, FDDObj)
    if FDDObj.fit_args.fddType == "FDD"; M = FDD(FDDObj.roi_xdat, x(1), x(2));
    elseif FDDObj.fit_args.fddType == "FDDG"; M = FDDG(FDDObj.roi_xdat, x(1), x(2), x(3), x(4));
    elseif FDDObj.fit_args.fddType == "FDDGpL"; M = FDDGpL(FDDObj.roi_xdat, x(1), x(2), x(3), x(4), x(5), x(6));
    elseif FDDObj.fit_args.fddType == "FDDGsL"; M = FDDGsL(FDDObj.roi_xdat, x(1), x(2), x(3), x(4), x(5), x(6));
    end
    M = x(7).*M;
    M(isnan(M)) = 0; if size(M, 2) > 1; M = M'; end
end

%% DEFINING THE FINAL FUNCTION THAT DESCRIBES THE FDD CURVE TO BE FITTED
function M = full_fdd_function(x, xdat, FDDObj)
    if FDDObj.fit_args.fddType == "FDD"; M = FDD(xdat, x(1), x(2));
    elseif FDDObj.fit_args.fddType == "FDDG"; M = FDDG(xdat, x(1), x(2), x(3), x(4));
    elseif FDDObj.fit_args.fddType == "FDDGpL"; M = FDDGpL(xdat, x(1), x(2), x(3), x(4), x(5), x(6));
    elseif FDDObj.fit_args.fddType == "FDDGsL"; M = FDDGsL(xdat, x(1), x(2), x(3), x(4), x(5), x(6));
    end
    M = x(7).*M;
    M(isnan(M)) = 0; if size(M, 2) > 1; M = M'; end
end