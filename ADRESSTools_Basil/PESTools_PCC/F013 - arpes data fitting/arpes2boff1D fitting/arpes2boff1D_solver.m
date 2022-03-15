function fitStr = arpes2boff1D_solver(arpesStr, cTYPE, FUNC, iparams, bTYPE, ibgrnd, solve_type)
% fitStr = arpes2boff1D_solver(arpesStr, cTYPE, FUNC, iparams, bTYPE, ibgrnd, solve_type)
%   Function that is used to fit N-parabolic dispersions to an EDC cut
%   through ARPES data based on a set of model functions (FUNC) that
%   contain the band offset vs subband energy. This allows the band offset
%   to be determined by extracting the best fit to the subbands,
%   constraining their energies to what is found in a theoretical model.
%   
%   IN:
%   -   arpesStr:       MATLAB data-structure that contains the initial ARPES data.
%   -   cTYPE:          1xN vector of the type of curve to use for fitting. Default: "sGLA" ("pGLA", "DS")
%   -   FUNC:           1xN cell of functions that gives the band offset vs subband energy
%   -   iparams:      	3 cells {x0}{lb}{ub} with 1x6 array: the model fit parameters [FDEF,FDT,FDW,BOFF,INT,FWHM]
%   -   bTYPE:          string of the type of background to use for fitting. Default: "Poly" ("none", "Shir", "LinShir")
%   -   ibgrnd:         3 cells {x0}{lb}{ub} with 1x5 vectors: the background parameters: [LHS,RHS,ORD,LAM,DEL,BGR]
%   -   solve_type:     String of either "fminunc", "fmincon", "lsqnonlin", "simulannealbnd" that sets the optimisation routine to use
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information

%% Default parameters
% - Default input parameters to use
if nargin < 7; solve_type = "fmincon"; end
if nargin < 6
    ibgrnd{1} = [min(arpesStr.xdat(:))+0.2, max(arpesStr.xdat(:))-0.2, 1, 0.5, 0, 0];
    ibgrnd{2} = abs(0.0.*ibgrnd{1}); ibgrnd{2}(:,6) = ibgrnd{1}(:,6) - 0.05;
    ibgrnd{3} = abs(0.0.*ibgrnd{1}); ibgrnd{2}(:,6) = ibgrnd{1}(:,6) + 0.05;
end
if nargin < 5; bTYPE = "Poly"; end
if isempty(solve_type); solve_type = "fmincon"; end
if isempty(ibgrnd)
    ibgrnd{1} = [min(arpesStr.xdat(:))+0.2, max(arpesStr.xdat(:))-0.2, 1, 0.5, 0, 0];
    ibgrnd{2} = abs(0.0.*ibgrnd{1}); ibgrnd{2}(:,6) = ibgrnd{1}(:,6) - 0.05;
    ibgrnd{3} = abs(0.0.*ibgrnd{1}); ibgrnd{2}(:,6) = ibgrnd{1}(:,6) + 0.05;
end
if isempty(bTYPE); bTYPE = "Poly"; end

%% - 1 - CONSTRAINING THE FIT PARAMETERS TO PHYSICAL VALUES
% -- Total number of states to be fitted
nSTATES	= length(cTYPE);
% -- Defining the initial conditions and limits of the optimisation method
x0 = [iparams{1}(:); ibgrnd{1}(:)];
lb = [iparams{2}(:); ibgrnd{2}(:)];
ub = [iparams{3}(:); ibgrnd{3}(:)];
% -- Applying constraints on the parameters, based on the data
% --- parameters
if x0(2) < 0; x0(2) = 0; end        % if FDT < 0, then make it 0
if x0(3) < 0; x0(3) = 0; end        % if FDW < 0, then make it 0
if x0(5) < 0; x0(5) = 0; end        % if MR < 0, then make it 0
if x0(5) > 1; x0(5) = 1; end        % if MR > 0, then make it 1
% --- lower bounds
if lb(2) < 0; lb(2) = 0; end        % if FDT < 0, then make it 0
if lb(3) < 0; lb(3) = 0; end        % if FDW < 0, then make it 0
if lb(5) < 0; lb(5) = 0; end        % if MR < 0, then make it 0
if lb(5) > 1; lb(5) = 1; end        % if MR > 0, then make it 1
% --- upper bounds
if ub(2) < 0; ub(2) = 0; end        % if FDT < 0, then make it 0
if ub(3) < 0; ub(3) = 0; end        % if FDW < 0, then make it 0
if ub(5) < 0; ub(5) = 0; end        % if MR < 0, then make it 0
if ub(5) > 1; ub(5) = 1; end        % if MR > 0, then make it 1
for i = 1:nSTATES
    % --- parameters
    if x0(5+i) < 0; x0(5+i) = 0; end	% if INT < 0, then make it 0
    if x0(8+i) < 0; x0(8+i) = 0; end    % if FWHM < 0, then make it 0
    % --- lower bounds
    if lb(5+i) < 0; lb(5+i) = 0; end   	% if INT < 0, then make it 0
    if lb(8+i) < 0; lb(8+i) = 0; end    % if FWHM < 0, then make it 0
    % --- upper bounds
    if ub(5+i) < 0; ub(5+i) = 0; end  	% if INT < 0, then make it 0
    if ub(8+i) < 0; ub(8+i) = 0; end    % if FWHM < 0, then make it 0
end
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
% -- Determining the degrees of freedom for the fitting model
val	= lb + (ub-lb); DoF = 0;
for i = 1:length(val); if val(i) ~= lb(i); DoF = DoF + 1; end; end

%% - 2 - DEFINING THE EDC DATA AND FITTING ARGUMENTS
% -- Extracting the data to be fitted
[roi_xdat, roi_ydat, roi_bgrnd] = PESBackground(arpesStr.xdat, arpesStr.ydat, bTYPE, x0(end-5), x0(end-4), x0(end-3), x0(end-2), x0(end-1), x0(end));
% -- Defining a structure that stores all relevant model and data variables
ARPESObj                = arpesStr;
ARPESObj.nSTATES        = nSTATES;
ARPESObj.roi_xdat       = roi_xdat;
ARPESObj.roi_ydat       = roi_ydat;
ARPESObj.roi_bgrnd      = roi_bgrnd;
ARPESObj.roi_ydat_sbtr	= ARPESObj.roi_ydat - ARPESObj.roi_bgrnd;
% -- Appending the input arguments to the global variable
ARPESObj.fit_args.solve_type    = solve_type;
ARPESObj.fit_args.FUNC          = FUNC;
ARPESObj.fit_args.cTYPE         = cTYPE;
ARPESObj.fit_args.iparams       = iparams;
ARPESObj.fit_args.bTYPE         = bTYPE;
ARPESObj.fit_args.ibgrnd        = ibgrnd;

%% - 3 - RUNNING THE SOLVER
% -- Local solver properties
MaxFunEvals = 1e5;
MaxIter     = 5e3;
FuncTol     = 1e-8; % 1e-7;
StepTol     = 1e-8; % 1e-7;
OptTol      = 1e-8; % 1e-7;
ConTol      = 1e-8; % 1e-7;
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
    [params,rnorm,resid,exitflag,output,lambda,jacobian]  = lsqcurvefit(@(x,roi_xdat) full_pes_function(x,roi_xdat,ARPESObj), x0, roi_xdat, roi_ydat, lb, ub, options);  
%% - 3.2 - LOCAL SOLVER: UNBOUNDED LEAST SQUARES NONLINEAR REGRESSION FIT
elseif solve_type == "nlinfit"
    % -- Defining the optimisation options for the simulated annealing method
    options = statset(...
        'MaxFunEvals', MaxFunEvals,...
        'MaxIter', MaxIter,...
        'TolFun', FuncTol);
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(roi_xdat, roi_ydat, @(x,roi_xdat) full_pes_function(x,roi_xdat,ARPESObj), x0, options);
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
    [params,fval,exitflag,output,grad,hessian] = fminunc(@(x) minimize_function(x,ARPESObj), x0, options); 
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
    [params,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) minimize_function(x,ARPESObj), x0, [], [], [], [], lb, ub, [], options);  
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
    [params,rnorm,resid,exitflag,output,lambda,jacobian] = lsqnonlin(@(x) minimize_function(x,ARPESObj), x0, lb, ub, options);    
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
    [params,fval,exitflag,output] = simulannealbnd(@(x) minimize_function(x,ARPESObj), x0, lb, ub, options);
end

%% - 4 - STORING THE FITS AS A MATLAB STRUCTURE
fitStr              = struct();
%% 4.1 - Storing the initial arguments
fitStr.solve_type	= solve_type;
fitStr.FUNC         = FUNC;
fitStr.cTYPE        = cTYPE;
fitStr.iparams      = iparams;
fitStr.bTYPE        = bTYPE;
fitStr.ibgrnd       = ibgrnd;
%% 4.2 - Storing the XPS data, background, model and minimisation variables
% -- Extracting the original data 
fitStr.xdat         = arpesStr.xdat;
fitStr.ydat          = arpesStr.ydat;
% -- Extracting the X-DOMAIN and DATA
[fitStr.X, fitStr.D] = fit_data(params, ARPESObj);
% -- Extracting the MODEL
fitStr.M        = fit_model(params, ARPESObj);
% -- Extracting the BACKGROUND
fitStr.B        = fit_background(params, ARPESObj);
% -- Extracting the DATA - BACKGROUND
fitStr.DB       = fitStr.D - fitStr.B;
fitStr.MB       = fitStr.M + fitStr.B;
% -- Extracting the RESIDUALS (R = D - (M + B))
fitStr.R        = fitStr.D - (fitStr.M + fitStr.B);
% -- Extracting the CHI-SQUARED value to be minimised
fitStr.CHISQ    = sum(fitStr.R.^2 ./ abs(fitStr.M + fitStr.B));
fitStr.MINFUN   = minimize_function(params, ARPESObj);
% -- Append the degrees of freedom
fitStr.DoF   	= DoF;
%% 4.3 - Storing the final fit variables for the PES background
fitStr.params	= params;
fitStr.nSTATES	= length(cTYPE);
fitStr.cTYPE  	= cTYPE;
% - Extracting parabolic parameters
fitStr.FDEF    = params(1);        % scalar of the FDD Fermi-Level position.
fitStr.FDT     = params(2);        % scalar of the FDD temperature.
fitStr.FDW     = params(3);        % scalar of the FDD Gaussian width after convolution.
fitStr.BOFF    = params(4);        % scalar of the band offset
fitStr.MR      = params(5);        % scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
fitStr.INT     = params(6:6+fitStr.nSTATES-1);                  	% scalar of the peak intensity of PE curve.
fitStr.FWHM    = params(6+fitStr.nSTATES:6+2*fitStr.nSTATES-1);       % scalar of the FWHM of the PE curve.
fitStr.LHS     = params(end-5);
fitStr.RHS     = params(end-4);
fitStr.ORD     = params(end-3);
fitStr.LAM     = params(end-2);
fitStr.DEL     = params(end-1);
fitStr.BGR     = params(end);
%% 4.4 - Storing the final fit variables of each PES curve component
fitStr.XX      	= linspace(min(fitStr.X(:)), max(fitStr.X(:)), 1e3)';
fitStr.YY      	= zeros(size(fitStr.XX));
for i = 1:fitStr.nSTATES
    BE(i)   = FUNC{i}(fitStr.BOFF);
    % --- Best fit components on new domain
    fitStr.cYY(:,i) = PESCurve_FDD(fitStr.XX, fitStr.cTYPE(i),...
        BE(i), fitStr.INT(i), fitStr.FWHM(i), fitStr.MR, 0, 0, 0, 0,...
        fitStr.FDEF, fitStr.FDT, fitStr.FDW);
    fitStr.YY = fitStr.YY + fitStr.cYY(:,i);
    % --- Storing each curve parameter
    fitStr.BE(i)    = BE(i);
    % --- Storing the curve area information
    fitStr.AREA(i)	= trapz(fitStr.XX, fitStr.cYY(:,i));
end
% --- Storing the normalised area contribution for quantification
fitStr.AREA0        = fitStr.AREA ./ trapz(fitStr.XX, fitStr.YY);

end

%% DEFINING THE RESIDUAL MINIMISATION FUNCTION TO BE MINIMISED
function MINFUN = minimize_function(x, ARPESObj)
    % - 1 - Extracting the DATA
    [~, D] = fit_data(x, ARPESObj);
    % - 2 - Extracting the MODEL
    M = fit_model(x, ARPESObj);
    % - 3 - Extracting the BACKGROUND
    B = fit_background(x, ARPESObj);
    % - 4 - Extracting the RESIDUALS (residuals = data - (model + background))
    R   = D - (M + B);
    % - 5 - Extracting the CHI-SQUARED value to be minimised
    MINFUN = sum(R.^2 ./ (M + B));
end

%% DEFINING THE FUNCTION THAT EXTRACTS THE XPS DATA TO BE FITTED
function [X, D] = fit_data(x, ARPESObj)
    [X, D, ~] = PESBackground(ARPESObj.xdat, ARPESObj.ydat, ARPESObj.fit_args.bTYPE, x(end-5), x(end-4), x(end-3), x(end-2), x(end-1), x(end));
    D(isnan(D)) = 0;
    if size(D, 2) > 1; D = D'; end
end

%% DEFINING THE FUNCTION THAT DETERMINES THE TOTAL PES BACKGROUND
function B = fit_background(x, ARPESObj)
    [~, ~, B] = PESBackground(ARPESObj.xdat, ARPESObj.ydat, ARPESObj.fit_args.bTYPE, x(end-5), x(end-4), x(end-3), x(end-2), x(end-1), x(end));
    B(isnan(B)) = 0;
    if size(B, 2) > 1; B = B'; end
end

%% DEFINING THE FUNCTION THAT DETERMINES THE TOTAL PES CURVE FIT
function M = fit_model(x, ARPESObj)
    % -- Extracting
    FUNC  	= ARPESObj.fit_args.FUNC;
    cTYPE  	= ARPESObj.fit_args.cTYPE;
    % -- Extracting the input variables for the fit
    FDEF    = x(1);        % scalar of the FDD Fermi-Level position.
    FDT     = x(2);        % scalar of the FDD temperature.
    FDW     = x(3);        % scalar of the FDD Gaussian width after convolution.
    BOFF    = x(4);        % scalar of the band offset
    MR      = x(5);        % scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
    INT     = x(6:6+ARPESObj.nSTATES-1);                     % scalar of the peak intensity of PE curve.
    FWHM    = x(6+ARPESObj.nSTATES:6+2*ARPESObj.nSTATES-1);	 % scalar of the FWHM of the PE curve.
    % - 1 - Extract the curve component for each state
    M     = zeros(size(ARPESObj.roi_xdat));      	% XPS model data
    for i = 1:ARPESObj.nSTATES
        BE = FUNC{i}(BOFF);
        cYY(:,i) = PESCurve_FDD(ARPESObj.roi_xdat, cTYPE(i), BE, INT(i), FWHM(i), MR, 0, 0, 0, 0, FDEF, FDT, FDW);
        M = M + cYY(:,i);
    end
    M(isnan(M)) = 0;
    if size(M, 2) > 1; M = M'; end
end

%% DEFINING THE FINAL FUNCTION THAT DESCRIBES THE PES CURVE TO BE FITTED
function MB = full_pes_function(x, xdat, ARPESObj)
    % - 1 - Extract the curve component for each state
    FUNC  	= ARPESObj.fit_args.FUNC;
    cTYPE  	= ARPESObj.fit_args.cTYPE;
    % -- Extracting the input variables for the fit
    FDEF    = x(1);        % scalar of the FDD Fermi-Level position.
    FDT     = x(2);        % scalar of the FDD temperature.
    FDW     = x(3);        % scalar of the FDD Gaussian width after convolution.
    BOFF    = x(4);        % scalar of the band offset
    MR      = x(5);        % scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
    INT     = x(6:6+ARPESObj.nSTATES-1);                     % scalar of the peak intensity of PE curve.
    FWHM    = x(6+ARPESObj.nSTATES:6+2*ARPESObj.nSTATES-1);	 % scalar of the FWHM of the PE curve.
    % - 1 - Extract the curve component for each state
    M     = zeros(size(ARPESObj.roi_xdat));      	% XPS model data
    for i = 1:ARPESObj.nSTATES
        BE = FUNC{i}(BOFF);
        cYY(:,i) = PESCurve_FDD(xdat, cTYPE(i), BE, INT(i), FWHM(i), MR, 0, 0, 0, 0, FDEF, FDT, FDW);
        M = M + cYY(:,i);
    end
    M(isnan(M)) = 0;
    if size(M, 2) > 1; M = M'; end
    % - 3 - Determine the background to be used
    [~, ~, B] = PESBackground(xdat, ARPESObj.roi_ydat, ARPESObj.fit_args.bTYPE, x(end-5), x(end-4), x(end-3), x(end-2), x(end-1), x(end));
    % - 4 - Final output is the sum of the model and background
    MB = M + B;
end