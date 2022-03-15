function fitStr = kz2model_solver(xdat, ydat, cTYPE, iparams, solve_type)
% fitStr = kz2model_solver(xdat, ydat, cTYPE, iparams, solve_type)
%   Function that runs a a Global Optimisation Algorithm based on Simulated
%   Annealing to determine the best curve fit to the input XPS data and
%   model (initial conditions and bounds). The temperature of the Simulated
%   Annealing is set to walk a large span of the parameter space to ensure
%   the global minimum is obtained. The minimisation function is based on
%   minimising the standard deviation of the squared residuals. Although
%   the "simulannealbnd()" algorithm is the best due to its global
%   convergence, the option is here to also use "lsqcurvefit()" and
%   "lsqnonlin()" for fast convergence to the XPS curve fits. Use this
%   function to run the convergence algorithm only one. 
%   NOTE: use 'xps_solver_n_runs()' to run the algorithm N times, allowing
%   for quantification of fit parameter uncertainties.
%   
%   IN:
%   -   xdat:           1xN vector that contains the domain of the line profile.
%   -   ydat:           1xN vector that contains the intensity data of the line profile.
%   -   cTYPE:          string of the type of curves to use; "sGLA" ("pGLA", "sGL", "pGL", "G", "L", "DS").
%   -   iparams:        3 cells {x0}{lb}{ub} with Nx8 arrays: the n'th peak parameters [KZ0,INT,MFP,DZ,MR,SKZ0,SINT,ASY,BGR]
%   -   solve_type:     string of the optimisation algorithm to use; "lsqcurvefit" (very fast, good), "lsqnonlin" (fast, not good), "simulannealbnd" (slow, best). 
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information

%% Default parameters
% - Default input parameters to use
if nargin < 5; solve_type = "lsqcurvefit"; end

%% - 1 - CONSTRAINING THE FIT PARAMETERS TO PHYSICAL VALUES
% -- Defining the initial conditions and limits of the optimisation method
x0 = iparams{1}(:);
lb = iparams{2}(:);
ub = iparams{3}(:);
% -- Defining the initial conditions and limits of the optimisation method
% --- parameters
if x0(2) < 0; x0(2) = 0; end    % if INT < 0, then make it 0
if x0(3) < 0; x0(3) = 0; end    % if MFP FWHM < 0, then make it 0
if x0(4) < 0; x0(4) = 0; end    % if DZ FWHM < 0, then make it 0
if x0(5) < 0; x0(5) = 0; end    % if MR < 0, then make it 0
if x0(5) > 1; x0(5) = 1; end    % if MR > 1, then make it 1
if x0(7) < 0; x0(7) = 0; end    % if SINT < 0, then make it 0
if x0(8) < 0; x0(8) = 0; end    % if ASY < 0, then make it 0
% --- lower bounds
if lb(2) < 0; lb(2) = 0; end    % if INT < 0, then make it 0
if lb(3) < 0; lb(3) = 0; end    % if MFP FWHM < 0, then make it 0
if lb(4) < 0; lb(4) = 0; end    % if DZ FWHM < 0, then make it 0
if lb(5) < 0; lb(5) = 0; end    % if MR < 0, then make it 0
if lb(5) > 1; lb(5) = 1; end    % if MR > 1, then make it 1
if lb(7) < 0; lb(7) = 0; end    % if SINT < 0, then make it 0
if lb(8) < 0; lb(8) = 0; end    % if ASY < 0, then make it 0
% --- lower bounds
if ub(2) < 0; ub(2) = 0; end    % if INT < 0, then make it 0
if ub(3) < 0; ub(3) = 0; end    % if MFP FWHM < 0, then make it 0
if ub(4) < 0; ub(4) = 0; end    % if DZ FWHM < 0, then make it 0
if ub(5) < 0; ub(5) = 0; end    % if MR < 0, then make it 0
if ub(5) > 1; ub(5) = 1; end    % if MR > 1, then make it 1
if ub(7) < 0; ub(7) = 0; end    % if SINT < 0, then make it 0
if ub(8) < 0; ub(8) = 0; end    % if ASY < 0, then make it 0
% -- Determining the degrees of freedom for the fitting model
val	= lb + (ub-lb); DoF = 0;
for i = 1:length(val); if val(i) ~= lb(i); DoF = DoF + 1; end; end

%% - 2 - DEFINING THE XPS DATA AND FITTING ARGUMENTS
% -- Defining a structure that stores all relevant model and data variables
DataObj        	= struct();
DataObj.xdat  	= xdat;
DataObj.ydat	= ydat;
% -- Appending the input arguments to the global variable
DataObj.fit_args.solve_type     = solve_type;
DataObj.fit_args.cTYPE          = cTYPE;
DataObj.fit_args.iparams        = iparams;

%% - 3 - RUNNING THE SOLVER
% -- Local solver properties
MaxFunEvals  = 1e5;
MaxIter     = 5e3;
FuncTol     = 1e-7; % 1e-7;
StepTol     = 1e-7; % 1e-7;
OptTol      = 1e-7; % 1e-7;
ConTol      = 1e-7; % 1e-7;
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
        'MaxIterations', MaxIter);
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,rnorm,resid,exitflag,output,lambda,jacobian]  = lsqcurvefit(@(x,xdat) full_fit_function(x,xdat,DataObj), x0, xdat, ydat, lb, ub, options);  
%% - 3.2 - LOCAL SOLVER: UNBOUNDED LEAST SQUARES NONLINEAR REGRESSION FIT
elseif solve_type == "nlinfit"
    % -- Defining the optimisation options for the simulated annealing method
    options = statset(...
        'MaxFunEvals', MaxFunEvals,...
        'MaxIter', MaxIter,...
        'TolFun', FuncTol);
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(roi_xdat, roi_int, @(x,xdat) full_fit_function(x,xdat,DataObj), x0, options);
%% - 3.3 - LOCAL SOLVER: FIND MINIMUM OF UNBOUNDED MULTIVARIATE FUNCTION
elseif solve_type == "fminunc"
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions('fminunc',...
        'MaxFunctionEvaluations', MaxFunEvals,...
        'FunctionTolerance', FuncTol,...
        'StepTolerance', StepTol,...
        'OptimalityTolerance', OptTol,...
        'MaxIterations', MaxIter);
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,fval,exitflag,output,grad,hessian] = fminunc(@(x) minimize_function(x,DataObj), x0, options); 
%% - 3.4 - LOCAL SOLVER: FIND MINIMUM OF BOUNDED MULTIVARIATE FUNCTION
elseif solve_type == "fmincon"
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions('fmincon',...
        'MaxFunctionEvaluations', MaxFunEvals,...
        'FunctionTolerance', FuncTol,...
        'StepTolerance', StepTol,...
        'OptimalityTolerance', OptTol,...
        'MaxIterations', MaxIter,...
        'ConstraintTolerance', ConTol);
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) minimize_function(x,DataObj), x0, [], [], [], [], lb, ub, [], options);  
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
    [params,rnorm,resid,exitflag,output,lambda,jacobian] = lsqnonlin(@(x) minimize_function(x,DataObj), x0, lb, ub, options);    
%% - 3.6 - GLOBAL SOLVER: BOUNDED SIMULATED ANNEALING TO MINIMISE RESIDUALS
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
    [params,fval,exitflag,output] = simulannealbnd(@(x) minimize_function(x,DataObj), x0, lb, ub, options);
end

%% - 4 - STORING THE FITS AS A MATLAB STRUCTURE
fitStr              = struct();
%% 4.1 - Storing the initial arguments
fitStr.solve_type	= solve_type;
fitStr.cTYPE        = cTYPE;
fitStr.iparams      = iparams;
%% 4.2 - Storing the XPS data, background, model and minimisation variables
% -- Extracting the original data 
fitStr.xdat     = xdat;
fitStr.ydat     = ydat;
% -- Extracting the X-DOMAIN and DATA
[fitStr.X, fitStr.D] = fit_data(params, DataObj);
% -- Extracting the MODEL
[fitStr.M, fitStr.M0, fitStr.M1] = fit_model(params, DataObj);
% -- Extracting the BACKGROUND
fitStr.B        = fit_background(params, DataObj);
% -- Extracting the DATA - BACKGROUND
fitStr.DB       = fitStr.D - fitStr.B;
fitStr.MB       = fitStr.M + fitStr.B;
% -- Extracting the RESIDUALS (R = D - (M + B))
fitStr.R        = fitStr.D - (fitStr.M + fitStr.B);
% -- Extracting the CHI-SQUARED value to be minimised
fitStr.CHISQ    = sum(fitStr.R.^2 ./ abs(fitStr.M + fitStr.B));
fitStr.MINFUN   = minimize_function(params, DataObj);
% -- Append the degrees of freedom
fitStr.DoF   	= DoF;
%% 4.3 - Storing the final fit variables of each PES curve component
fitStr.XX    	= linspace(min(fitStr.X(:)), max(fitStr.X(:)), 1e3)';
fitStr.YY0    	= PESCurve(fitStr.XX, fitStr.cTYPE, params(1), params(2), 1./params(3), params(5), params(6), params(7), 0.0, params(8));
fitStr.YY1    	= PESCurve(fitStr.XX, fitStr.cTYPE, mean([params(1),params(1)+params(6)]), params(2), 1./params(4), params(5), 0.0, 0.0, 0.0, params(8));
fitStr.YY1      = (fitStr.YY1 ./ max(fitStr.YY1(:)))*max(fitStr.YY0(:));
fitStr.YY    	= full_fit_function(params, fitStr.XX, DataObj);
fitStr.params 	= params;
fitStr.KZ0      = params(1);
fitStr.INT      = params(2);
fitStr.MFP      = params(3);
fitStr.DZ       = params(4);
fitStr.MR       = params(5);
fitStr.SKZ0 	= params(6);
fitStr.SINT 	= params(7);
fitStr.ASY      = params(8);
fitStr.BGR      = params(9);

end

%% DEFINING THE RESIDUAL MINIMISATION FUNCTION TO BE MINIMISED
function MINFUN = minimize_function(x, DataObj)
    % - 1 - Extracting the DATA
    [~, D]  = fit_data(x, DataObj);
    % - 2 - Extracting the MODEL
    [M, ~, ~]  = fit_model(x, DataObj);
    % - 3 - Extracting the BACKGROUND
    B       = fit_background(x, DataObj);
    % - 4 - Extracting the RESIDUALS (residuals = data - (model + background))
    R       = D - (M + B);
    % - 5 - Extracting the CHI-SQUARED value to be minimised
    MINFUN  = sum(R.^2 ./ (M + B));
end

%% DEFINING THE FUNCTION THAT EXTRACTS THE DATA TO BE FITTED
function [X, D] = fit_data(x, DataObj)
    X   = DataObj.xdat;
    D   = DataObj.ydat;
    D(isnan(D)) = 0;
    if size(D, 2) > 1; D = D'; end
end

%% DEFINING THE FUNCTION THAT DETERMINES THE TOTAL PES BACKGROUND
function B = fit_background(x, DataObj)
    B           = DataObj.xdat.*0 + x(9);
    B(isnan(B)) = 0;
    if size(B, 2) > 1; B = B'; end
end

%% DEFINING THE FUNCTION THAT DETERMINES THE TOTAL PES CURVE FIT
function [M, M0, M1] = fit_model(x, DataObj)
    % - 1 - Extract the curve with only MFP broadening
    M0  = PESCurve(DataObj.xdat, DataObj.fit_args.cTYPE, x(1), x(2), 1./x(3), x(5), x(6), x(7), 0.0, x(8));
    M0(isnan(M0)) = 0;
    if size(M0, 2) > 1; M0 = M0'; end
    % - 2 - Extract the curve with only kz broadening
    M1  = PESCurve(DataObj.xdat, DataObj.fit_args.cTYPE, x(1), x(2), 1./x(4), x(5), x(6), x(7), 0.0, x(8));
    M1(isnan(M1)) = 0;
    if size(M1, 2) > 1; M1 = M1'; end
    % - 3 - Extract the curve with MFP + KZ broadening
    M  = PESCurve(DataObj.xdat, DataObj.fit_args.cTYPE, x(1), x(2), (1./x(3)+1./x(4)), x(5), x(6), x(7), 0.0, x(8));
    M(isnan(M)) = 0;
    if size(M, 2) > 1; M = M'; end
end

%% DEFINING THE FINAL FUNCTION THAT DESCRIBES THE PES CURVE TO BE FITTED
function MB = full_fit_function(x, xdat, DataObj)
    % - 1 - Extract the curve with MFP + KZ broadening
    M  = PESCurve(xdat, DataObj.fit_args.cTYPE, x(1), x(2), (1./x(3)+1./x(4)), x(5), x(6), x(7), 0.0, x(8));
    M(isnan(M)) = 0;
    if size(M, 2) > 1; M = M'; end
    % - 2 - Extracting the BACKGROUND    
    B = xdat.*0 + x(9);
    B(isnan(B)) = 0;
    if size(B, 2) > 1; B = B'; end
    % - 3 - Extracting the FULL MODEL CURVE
    MB      = M + B;
end