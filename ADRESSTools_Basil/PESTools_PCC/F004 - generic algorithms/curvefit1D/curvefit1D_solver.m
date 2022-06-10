function fitStr = curvefit1D_solver(xdat, ydat, curve_type, init_params, bgrnd_type, init_bgrnd, solve_type)
% fitStr = curvefit1D_solver(xdat, ydat, curve_type, init_params, bgrnd_type, init_bgrnd, solve_type)
%   Function that runs an Optimisation Algorithm to determine the best fit
%   for N-curves of a 1D data set defined by 'xdat' (domain) and 'ydat'
%   (range). If 'curve_type' is empty, an asymmetric pseudo-Voigt curve is
%   used. 'init_params' is the cell array that defines the initial guess
%   and lower-, upper-bound constraints of the fit parameters.
%   
%   IN:
%   -   xdat:           [nX x 1] array of the x-axis data
%   -   ydat:           [nY x 1] array of the y-axis data
%   -   curve_type:  	string of the type of curve to be used for fitting. Default: "sGLA" ("pGLA", "sGL", "pGL", "G", "L", "DS")
%   -   init_params:   	{3 x 1} cell array of {x0}{lb}{ub}, where each cell is an [N x 8] matrix of N'th peak parameters: [XLOC,YINT,FWHM,MR,ASY,LSX,LSY,LSW]
%   -   [bgrnd_type]:  	string of the type of background to use for fitting. Default: "none" ("Poly", "Shir", "LinShir")
%   -   [init_bgrnd]:   {3 x 1} cell array of {x0}{lb}{ub}, where each cell is an [1 x 5] matrix of the background parameters: [LHS,RHS,ORD,LAM,DEL,BGR]
%   -   [solve_type]:   string of the optimisation algorithm to use; "lsqcurvefit" (very fast, good), "fminunc" (quick, good), "simulannealbnd" (slow, best). 
%                       Default: "lsqcurvefit" ("lsqcurvefit", "fminunc", "fmincon", "fminunc-std", "simulannealbnd", "simulannealbnd-std")
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all fit parameters / variables / information

%% Default input parameters & consistency checks
pp  = plot_props();
% -- Defining input parameters if undefined
if nargin < 7; solve_type = "lsqcurvefit"; end
if nargin < 6; init_bgrnd{1} = [min(xdat(:)), max(xdat(:)), 1, 0.5, 0, 0]; init_bgrnd{2} = init_bgrnd{1}; init_bgrnd{2}(6) = -0.10; init_bgrnd{3} = init_bgrnd{1}; init_bgrnd{3}(6) = 0.10; end
if nargin < 5; bgrnd_type = "none"; end
if isempty(solve_type); solve_type = "lsqcurvefit"; end
if isempty(init_bgrnd); init_bgrnd{1} = [min(xdat(:)), max(xdat(:)), 1, 0.5, 0, 0]; init_bgrnd{2} = init_bgrnd{1}; init_bgrnd{2}(6) = -0.10; init_bgrnd{3} = init_bgrnd{1}; init_bgrnd{3}(6) = 0.10; end
if isempty(bgrnd_type); bgrnd_type = "none"; end
if isempty(curve_type); curve_type = "sGLA"; end
% -- Consistency check and finding the total number of curves
if size(init_params{1}, 1) ~= size(init_params{2}, 1) || size(init_params{2}, 1) ~= size(init_params{3}, 1)
    error('The input parameter cell array is not a consistent size - check iparams input!');
end
% -- Extracting the total number of curves
ncurves = size(init_params{1}, 1);
% -- Validity check on input parameters
for i = 1:ncurves
    if size(init_params{1}, 2) == 4
        % --- Set the ASY, LSE, LSI and LSW components to zero
        init_params{1}(i,5) = 0; init_params{1}(i,6) = 0; init_params{1}(i,7) = 0; init_params{1}(i,8) = 0;
        init_params{2}(i,5) = 0; init_params{2}(i,6) = 0; init_params{2}(i,7) = 0; init_params{2}(i,8) = 0;
        init_params{3}(i,5) = 0; init_params{3}(i,6) = 0; init_params{3}(i,7) = 0; init_params{3}(i,8) = 0;
    elseif size(init_params{1}, 2) == 5
        % --- Set the LSE, LSI and LSW components to zero
        init_params{1}(i,6) = 0; init_params{1}(i,7) = 0; init_params{1}(i,8) = 0;
        init_params{2}(i,6) = 0; init_params{2}(i,7) = 0; init_params{2}(i,8) = 0;
        init_params{3}(i,6) = 0; init_params{3}(i,7) = 0; init_params{3}(i,8) = 0;
    elseif size(init_params{1}, 2) == 6
        % --- Set the LSI and LSW components to zero
        init_params{1}(i,7) = 0; init_params{1}(i,8) = 0;
        init_params{2}(i,7) = 0; init_params{2}(i,8) = 0;
        init_params{3}(i,7) = 0; init_params{3}(i,8) = 0;
    elseif size(init_params{1}, 2) == 7
        % --- Set the LSW components to zero
        init_params{1}(i,8) = 0;
        init_params{2}(i,8) = 0;
        init_params{3}(i,8) = 0;
    elseif size(init_params{1}, 2) > 8 || size(init_params{1}, 2) < 4 
        error('Not enough input arguments defined - check iparams input!');
    end
end

%% - 1 - CONSTRAINING THE FIT PARAMETERS TO PHYSICAL VALUES
% -- Defining the initial conditions and limits of the optimisation method
x0 = [init_params{1}(:); init_bgrnd{1}(:)];
lb = [init_params{2}(:); init_bgrnd{2}(:)];
ub = [init_params{3}(:); init_bgrnd{3}(:)];
% -- Defining the initial conditions and limits of the optimisation method
for i = 1:ncurves
    indx = i:ncurves:(length(x0)-6);
    X0 = x0(indx); LB = lb(indx); UB = ub(indx);
    % --- parameters
    if X0(2) < 0; x0(indx(2)) = 0; end    % if YINT < 0, then make it 0
    if X0(3) < 0; x0(indx(3)) = 0; end    % if FWHM < 0, then make it 0
    if X0(4) < 0; x0(indx(4)) = 0; end    % if MR < 0, then make it 0
    if X0(4) > 1; x0(indx(4)) = 1; end    % if MR > 1, then make it 1
    if X0(5) < 0; x0(indx(5)) = 0; end    % if ASY < 0, then make it 0
    if X0(7) < 0; x0(indx(7)) = 0; end    % if LSY < 0, then make it 0
    if X0(8) < 0; x0(indx(8)) = 0; end    % if LSW < 0, then make it 0
    % --- lower bounds
    if LB(2) < 0; lb(indx(2)) = 0; end    % if YINT < 0, then make it 0
    if LB(3) < 0; lb(indx(3)) = 0; end    % if FWHM < 0, then make it 0
    if LB(4) < 0; lb(indx(4)) = 0; end    % if MR < 0, then make it 0
    if LB(4) > 1; lb(indx(4)) = 1; end    % if MR > 1, then make it 1
    if LB(5) < 0; lb(indx(5)) = 0; end    % if ASY < 0, then make it 0
    if LB(7) < 0; lb(indx(7)) = 0; end    % if LSY < 0, then make it 0
    if LB(8) < 0; lb(indx(8)) = 0; end    % if LSW < 0, then make it 0
    % --- lower bounds
    if UB(2) < 0; ub(indx(2)) = 0; end    % if YINT < 0, then make it 0
    if UB(3) < 0; ub(indx(3)) = 0; end    % if FWHM < 0, then make it 0
    if UB(4) < 0; ub(indx(4)) = 0; end    % if MR < 0, then make it 0
    if UB(4) > 1; ub(indx(4)) = 1; end    % if MR > 1, then make it 1
    if UB(5) < 0; ub(indx(5)) = 0; end    % if ASY < 0, then make it 0
    if UB(7) < 0; ub(indx(7)) = 0; end    % if LSY < 0, then make it 0
    if UB(8) < 0; ub(indx(8)) = 0; end    % if LSW < 0, then make it 0
    % --- forcing the background polynomial order to be an integer
    x0(end-3) = round(ceil(x0(end-3)));
    lb(end-3) = round(ceil(lb(end-3)));
    ub(end-3) = round(ceil(ub(end-3)));
    % --- forcing the lambda parameters between 0->1 (linear vs shirley mixing ratio)
    if x0(end-2) < 0; x0(end-2) = 0; end	% if LAM < 0, then make it 0
    if lb(end-2) < 0; lb(end-2) = 0; end	% if LAM < 0, then make it 0
    if ub(end-2) < 0; ub(end-2) = 0; end	% if LAM < 0, then make it 0
    if x0(end-2) > 1; x0(end-2) = 1; end	% if LAM > 1, then make it 1
    if lb(end-2) > 1; lb(end-2) = 1; end	% if LAM > 1, then make it 1
    if ub(end-2) > 1; ub(end-2) = 1; end	% if LAM > 1, then make it 1    
end
% -- Determining the degrees of freedom for the fitting model
val	= lb + (ub-lb); DoF = 0;
for i = 1:length(val); if val(i) ~= lb(i); DoF = DoF + 1; end; end

%% - 2 - DEFINING THE ROI DATA AND FITTING ARGUMENTS
% -- Extracting the data to be fitted
[roi_xdat, roi_ydat, roi_bgrnd] = PESBackground(xdat, ydat, bgrnd_type, x0(end-5), x0(end-4), x0(end-3), x0(end-2), x0(end-1), x0(end));
% -- Defining a structure that stores all relevant model and data variables
FitObj                          = struct();
FitObj.xdat                     = xdat;
FitObj.ydat                     = ydat;
FitObj.roi_xdat                 = roi_xdat;
FitObj.roi_ydat                 = roi_ydat;
FitObj.roi_bgrnd                = roi_bgrnd;
% -- Appending the input arguments to the global variable
FitObj.fit_args.solve_type      = solve_type;
FitObj.fit_args.ncurves         = ncurves;
FitObj.fit_args.curve_type    	= curve_type;
FitObj.fit_args.init_params     = init_params;
FitObj.fit_args.bgrnd_type    	= bgrnd_type;
FitObj.fit_args.init_bgrnd      = init_bgrnd;

%% - 3 - RUNNING THE SOLVER
% -- Local solver properties
MaxFunEvals = 1e5;
MaxIter     = 5e3;
FuncTol     = 1e-10; % 1e-7;
StepTol     = 1e-10; % 1e-7;
OptTol      = 1e-10; % 1e-7;
ConTol      = 1e-10; % 1e-7;
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
        'MaxIterations', MaxIter); % 'Display', 'iter'
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,rnorm,resid,exitflag,output,lambda,jacobian]  = lsqcurvefit(@(x,roi_xdat) full_pes_function(x,roi_xdat,FitObj), x0, roi_xdat, roi_ydat, lb, ub, options);  

%% - 3.2 - LOCAL SOLVER: FIND MINIMUM OF UNBOUNDED MULTIVARIATE FUNCTION
elseif solve_type == "fminunc"
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions('fminunc',...
        'MaxFunctionEvaluations', MaxFunEvals,...
        'FunctionTolerance', FuncTol,...
        'StepTolerance', StepTol,...
        'OptimalityTolerance', OptTol,...
        'MaxIterations', MaxIter); % 'Display', 'iter'
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,fval,exitflag,output,grad,hessian] = fminunc(@(x) minimize_function(x,FitObj), x0, options); 

%% - 3.3 - LOCAL SOLVER: FIND MINIMUM OF BOUNDED MULTIVARIATE FUNCTION
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
    [params,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) minimize_function(x,FitObj), x0, [], [], [], [], lb, ub, [], options);  

%% - 3.4 - LOCAL SOLVER: FIND MINIMUM OF BOUNDED MULTIVARIATE FUNCTION (MINIMISE STD)
elseif solve_type == "fmincon-std"
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions('fmincon',...
        'MaxFunctionEvaluations', MaxFunEvals,...
        'FunctionTolerance', FuncTol,...
        'StepTolerance', StepTol,...
        'OptimalityTolerance', OptTol,...
        'MaxIterations', MaxIter,...
        'ConstraintTolerance', ConTol);
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) minimize_function_std(x,FitObj), x0, [], [], [], [], lb, ub, [], options);  

    
%% - 3.5 - GLOBAL SOLVER: BOUNDED SIMULATED ANNEALING TO MINIMISE RESIDUALS
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
    [params,fval,exitflag,output] = simulannealbnd(@(x) minimize_function(x,FitObj), x0, lb, ub, options);
    
%% - 3.6 - GLOBAL SOLVER: BOUNDED SIMULATED ANNEALING TO MINIMISE STD
elseif solve_type == "simulannealbnd-std"
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
    [params,fval,exitflag,output] = simulannealbnd(@(x) minimize_function_std(x,FitObj), x0, lb, ub, options);
end

%% - 4 - STORING THE FITS AS A MATLAB STRUCTURE
fitStr              = struct();
%% 4.1 - Storing the initial arguments
fitStr.solve_type	= solve_type;
fitStr.curve_type	= curve_type;
fitStr.init_params	= init_params;
fitStr.bgrnd_type	= bgrnd_type;
fitStr.init_bgrnd 	= init_bgrnd;
fitStr.ncurves      = ncurves;
%% 4.2 - Storing the XPS data, background, model and minimisation variables
% -- Extracting the original data 
fitStr.xdat     = xdat;
fitStr.ydat 	= ydat;
% -- Extracting the X-DOMAIN and DATA
[fitStr.X, fitStr.D] = fit_data(params, FitObj);
% -- Extracting the MODEL
fitStr.M        = fit_model(params, FitObj);
% -- Extracting the BACKGROUND
fitStr.B        = fit_background(params, FitObj);
% -- Extracting the DATA - BACKGROUND
fitStr.DB       = fitStr.D - fitStr.B;
fitStr.MB       = fitStr.M + fitStr.B;
% -- Extracting the RESIDUALS (R = M - (D + B))
fitStr.R        = fitStr.D - (fitStr.M + fitStr.B);
% -- Extracting the CHI-SQUARED value to be minimised
fitStr.CHISQ    = sum(fitStr.R.^2 ./ abs(fitStr.M + fitStr.B));
fitStr.MINFUN   = minimize_function(params, FitObj);
% -- Append the degrees of freedom
fitStr.DoF   	= DoF;
%% 4.3 - Storing the final fit variables for the PES background
fitStr.bgrnd_params	= params(end-5:end);
fitStr.LHS          = fitStr.bgrnd_params(1);
fitStr.RHS          = fitStr.bgrnd_params(2);
fitStr.ORD          = fitStr.bgrnd_params(3);
fitStr.LAM          = fitStr.bgrnd_params(4);
fitStr.DEL          = fitStr.bgrnd_params(5);
fitStr.BGR          = fitStr.bgrnd_params(6);
%% 4.4 - Storing the final fit variables of each PES curve component
fitStr.XX           = linspace(min(fitStr.X(:)), max(fitStr.X(:)), 1e3)';
fitStr.YY           = zeros(size(fitStr.XX));
for i = 1:fitStr.ncurves
    % --- Best fit components on new domain
    fitStr.curve_params(i,:) = params(i:fitStr.ncurves:(end-6));
   	fitStr.cYY(:,i) = PESCurve(fitStr.XX, fitStr.curve_type,...
        fitStr.curve_params(i,1), fitStr.curve_params(i,2), fitStr.curve_params(i,3),...
        fitStr.curve_params(i,4), fitStr.curve_params(i,6), fitStr.curve_params(i,7),...
        fitStr.curve_params(i,8), fitStr.curve_params(i,5));
    fitStr.YY = fitStr.YY + fitStr.cYY(:,i);
    % --- Storing each curve parameter
    fitStr.XLOC(i)  = fitStr.curve_params(i,1);
    fitStr.YINT(i)  = fitStr.curve_params(i,2);
    fitStr.FWHM(i)  = fitStr.curve_params(i,3);
    fitStr.MR(i)    = fitStr.curve_params(i,4);
    fitStr.ASY(i)   = fitStr.curve_params(i,5);
    fitStr.LSX(i)   = fitStr.curve_params(i,6);
    fitStr.LSY(i)   = fitStr.curve_params(i,7);
    fitStr.LSW(i)   = fitStr.curve_params(i,8);
    % --- Storing the curve area information
    fitStr.AREA(i)	= trapz(fitStr.XX, fitStr.cYY(:,i));
end
% --- Storing the normalised area contribution for quantification
fitStr.AREA0        = fitStr.AREA ./ trapz(fitStr.XX, fitStr.YY);

end

%% DEFINING THE RESIDUAL MINIMISATION FUNCTION TO XLOC MINIMISED
function MINFUN = minimize_function(x, FitObj)
    % - 1 - Extracting the DATA
    [~, D] = fit_data(x, FitObj);
    % - 2 - Extracting the MODEL
    M = fit_model(x, FitObj);
    % - 3 - Extracting the BACKGROUND
    B = fit_background(x, FitObj);
    % - 4 - Extracting the RESIDUALS (residuals = data - (model + background))
    R = D - (M + B);
    % - 5 - Extracting the CHI-SQUARED value to be minimised
    MINFUN = sum(R.^2 ./ abs(M + B));
end

%% DEFINING THE FUNCTION TO XLOC MINIMISED VIA STD
function MINFUN = minimize_function_std(x, FitObj)
    % - 1 - Extracting the DATA
    [~, D] = fit_data(x, FitObj);
    % - 2 - Extracting the MODEL
    M   = fit_model(x, FitObj);
    % - 3 - Extracting the BACKGROUND
    B   = fit_background(x, FitObj);
    % - 4 - Extracting the RESIDUALS (residuals = data - (model + background))
    R = D - (M + B);
    % - 5 - Extracting the CHI-SQUARED value to be minimised
    MINFUN = std(R.^2 ./ abs(M + B));
end

%% DEFINING THE FUNCTION THAT EXTRACTS THE XPS DATA TO XLOC FITTED
function [X, D] = fit_data(x, FitObj)
    [X, D, ~] = PESBackground(FitObj.xdat, FitObj.ydat, FitObj.fit_args.bgrnd_type, x(end-5), x(end-4), x(end-3), x(end-2), x(end-1), x(end));
    D(isnan(D)) = 0;
    if size(D, 2) > 1; D = D'; end
end

%% DEFINING THE FUNCTION THAT DETERMINES THE TOTAL PES BACKGROUND
function B = fit_background(x, FitObj)
    [~, ~, B] = PESBackground(FitObj.xdat, FitObj.ydat, FitObj.fit_args.bgrnd_type, x(end-5), x(end-4), x(end-3), x(end-2), x(end-1), x(end));
    B(isnan(B)) = 0;
    if size(B, 2) > 1; B = B'; end
end

%% DEFINING THE FUNCTION THAT DETERMINES THE TOTAL PES CURVE FIT
function M = fit_model(x, FitObj)
    % - 1 - Extract the curve component for each state
    comp_int = {};
    for i = 1:FitObj.fit_args.ncurves
        % -- Extracting the arguments for the component curve
        pes_args    = x(i:FitObj.fit_args.ncurves:(end-6));
        % -- Extracting the component intensities
        comp_int{i} = PESCurve(FitObj.roi_xdat, FitObj.fit_args.curve_type,...
            pes_args(1), pes_args(2), pes_args(3),...
            pes_args(4), pes_args(6), pes_args(7),...
            pes_args(8), pes_args(5));
    end
    % - 2 - Sum up all of the curve components to get the final model PES curve
    pes_int = zeros(size(FitObj.roi_xdat));
    for i = 1:FitObj.fit_args.ncurves
        pes_int = pes_int + comp_int{i};
    end
    M = pes_int;
    M(isnan(M)) = 0;
    if size(M, 2) > 1; M = M'; end
end

%% DEFINING THE FINAL FUNCTION THAT DESCRIXLOCS THE PES CURVE TO XLOC FITTED
function MB = full_pes_function(x, xdat, FitObj)
    % - 1 - Extract the curve component for each state
    comp_int = {};
    for i = 1:FitObj.fit_args.ncurves
        % -- Extracting the arguments for the component curve
        pes_args    = x(i:FitObj.fit_args.ncurves:(end-6));
        % -- Extracting the component intensities
        comp_int{i} = PESCurve(xdat, FitObj.fit_args.curve_type,...
            pes_args(1), pes_args(2), pes_args(3),...
            pes_args(4), pes_args(6), pes_args(7),...
            pes_args(8), pes_args(5));
    end
    % - 2 - Sum up all of the curve components to get the final model PES curve
    pes_int = zeros(size(xdat));
    for i = 1:FitObj.fit_args.ncurves
        pes_int = pes_int + comp_int{i};
    end
    M = pes_int;
    % - 3 - Determine the background to be used
    [~, ~, B] = PESBackground(xdat, FitObj.roi_ydat, FitObj.fit_args.bgrnd_type, x(end-5), x(end-4), x(end-3), x(end-2), x(end-1), x(end));
    % - 4 - Final output is the sum of the model and background
    MB = M + B;
end