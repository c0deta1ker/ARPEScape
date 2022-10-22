function fitStr = arpes2boff2D_solver(arpesStr, modelFunc, cTYPE, iparams, ibgrnd, solve_type)
% fitStr = arpes2boff2D_solver(arpesStr, modelFunc, cTYPE, iparams, ibgrnd, solve_type)
%   Function that is used to fit N-parabolic dispersions to ARPES data, by
%   defining the energy/momentum minimum and gaussian broadening, and
%   effective mass of each subband state. A Fermi-Dirac Distribution is
%   also used with a linear background, for states that are close to the
%   Fermi edge. Constrained so that the subband energies are given by a
%   functional form modelFunc.
%
%   IN:
%   -   arpesStr:       MATLAB data-structure that contains the initial ARPES data.
%   -   modelFunc:      1xN cell of functions that gives the band offset vs subband energy
%   -   cTYPE:          1xN vector of the type of curve to use for fitting. Default: "G2DA" ("G2D", "L2D")
%   -   iparams:      	3 cells {x0}{lb}{ub} with 1x6 array: the model fit parameters [INT,XLOC,YLOC,XFWHM,YFWHM,MSTAR]
%   -   ibgrnd:         3 cells {x0}{lb}{ub} with 1x6 vector of the background parameters: [FDEF,FDT,FDW,BGR,BIN,BCO]
%   -   solve_type:     String of either "fminunc", "fmincon", "lsqnonlin", "simulannealbnd" that sets the optimisation routine to use
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information

%% Default parameters & validity checks
% - Default input parameters to use
if nargin < 6; solve_type = "fmincon"; end
if isempty(solve_type); solve_type = "fmincon"; end

%% - 1 - CONSTRAINING THE FIT PARAMETERS TO PHYSICAL VALUES
% -- Extracting the total number of states to be fitted
n       = length(cTYPE);
% -- Defining the initial conditions and limits of the parameters
x0      = [iparams{1}(:); ibgrnd{1}(:)];
lb      = [iparams{2}(:); ibgrnd{2}(:)];
ub      = [iparams{3}(:); ibgrnd{3}(:)];    
% -- Applying constraints on the parameters, based on the data
for i = 1:n
    % --- parameters
    if x0(i+2*n-1) < 0; x0(i+2*n-1) = 0; end   	% if INT < 0, then make it 0
    if x0(i+3*n-1) < 0; x0(i+3*n-1) = 0; end    % if XFWHM < 0, then make it 0
    if x0(i+4*n-1) < 0; x0(i+4*n-1) = 0; end    % if YFWHM < 0, then make it 0
    % --- lower bounds
    if lb(i+2*n-1) < 0; lb(i+2*n-1) = 0; end  	% if INT < 0, then make it 0
    if lb(i+3*n-1) < 0; lb(i+3*n-1) = 0; end    % if XFWHM < 0, then make it 0
    if lb(i+4*n-1) < 0; lb(i+4*n-1) = 0; end    % if YFWHM < 0, then make it 0
    % --- upper bounds
    if ub(i+2*n-1) < 0; ub(i+2*n-1) = 0; end  	% if INT < 0, then make it 0
    if ub(i+3*n-1) < 0; ub(i+3*n-1) = 0; end    % if XFWHM < 0, then make it 0
    if ub(i+4*n-1) < 0; ub(i+4*n-1) = 0; end    % if YFWHM < 0, then make it 0
end 
% -- Determining the degrees of freedom for the fitting model
val	= lb + (ub-lb); DoF = 0;
for i = 1:length(val); if val(i) ~= lb(i); DoF = DoF + 1; end; end

%% - 2 - DEFINING THE ARPES DATA AND FITTING ARGUMENTS
% -- Defining global variable that stores all relevant model and data variables
ARPESObj            = struct();
ARPESObj.arpesStr	= arpesStr;
% -- Appending the input arguments to the global variable
ARPESObj.fit_args.solve_type    = solve_type;
ARPESObj.fit_args.cTYPE         = cTYPE;
ARPESObj.fit_args.modelFunc          = modelFunc;
ARPESObj.fit_args.iparams   	= iparams;
ARPESObj.fit_args.ibgrnd    	= ibgrnd;

%% - 3 - RUNNING THE SOLVER
% -- Local solver properties
MaxFunEvals = 1e5;
MaxIter     = 1e4;
FuncTol     = 1e-12;
StepTol     = 1e-12;
OptTol      = 1e-12;
ConTol      = 1e-12;
% -- Simulated annealing properties
TempFcn     = 'temperaturefast';
InitTemp    = 1e2;
ReannealInt = 1e2;
%% - 3.1 - LOCAL SOLVER: FIND MINIMUM OF UNBOUNDED MULTIVARIATE FUNCTION
if solve_type == "fminunc"
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions('fminunc',...
        'MaxFunctionEvaluations', MaxFunEvals,...
        'FunctionTolerance', FuncTol,...
        'StepTolerance', StepTol,...
        'OptimalityTolerance', OptTol,...
        'MaxIterations', MaxIter);
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,fval,exitflag,output,lambda,grad,hessian] = fminunc(@(x) minimize_function(x,ARPESObj), x0, options); 
%% - 3.2 - LOCAL SOLVER: FIND MINIMUM OF BOUNDED MULTIVARIATE FUNCTION
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
%% - 3.3 - LOCAL SOLVER: BOUNDED LEAST SQUARES FITTING METHOD FOR NON-LINEAR LEAST SQUARES PROBLEMS
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
%% - 3.4 - GLOBAL SOLVER: BOUNDED SIMULATED ANNEALING TO MINIMISE RESIDUALS
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
fitStr.modelFunc         = modelFunc;
fitStr.iparams      = iparams;
fitStr.ibgrnd       = ibgrnd;
%% 4.2 - Storing the ARPES data, background, model and minimisation variables
[ARPES, MODEL]      = extract_arpes_and_model(params, ARPESObj.arpesStr, cTYPE, modelFunc);
% -- Storing the domain and range
fitStr.kx           = MODEL.kx;
fitStr.eb           = MODEL.eb;
% -- x- and y-axis limits
fitStr.kx_lims      = [min(fitStr.kx(:)), max(fitStr.kx(:))];
fitStr.eb_lims      = [min(fitStr.eb(:)), max(fitStr.eb(:))];
% -- Extracting the  DATA
fitStr.D        = ARPES.data;
% -- Extracting the MODEL
fitStr.M        = MODEL.data;
% -- Extracting the RESIDUALS (R = D - M)
fitStr.R        = fitStr.D - fitStr.M;
fitStr.Rx    	= mean(fitStr.R, 1);      % R in 1D along x
fitStr.Ry       = mean(fitStr.R, 2);      % R in 1D along y
% -- Extracting the CHI-SQUARED value to be minimised
fitStr.CHISQ2D	= fitStr.R.^2 ./ abs(fitStr.M);
fitStr.CHISQx 	= mean(fitStr.CHISQ2D, 1);      % CHISQ in 1D along x
fitStr.CHISQy  	= mean(fitStr.CHISQ2D, 2);      % CHISQ in 1D along y
fitStr.CHISQ    = sum(sum(fitStr.CHISQ2D));
fitStr.MINFUN   = minimize_function(params, ARPESObj);
% -- Append the degrees of freedom
fitStr.DoF   	= DoF;
%% 4.3 - Storing the final fit variables for the ARPES model fit
fitStr.params	= params;
fitStr.nSTATES  = length(cTYPE);
fitStr.cTYPE  	= cTYPE;
% - Extracting parabolic parameters
fitStr.BOFF    = params(1);
fitStr.XLOC    = params(2:2+n-1);
fitStr.YLOC = []; for i = 1:n; fitStr.YLOC(i) = modelFunc{i}(fitStr.BOFF); end
fitStr.INT     = params(2+n:2+2*n-1); 
fitStr.XFWHM   = params(2+2*n:2+3*n-1);
fitStr.YFWHM   = params(2+3*n:2+4*n-1);
fitStr.MSTAR   = params(2+4*n:2+5*n-1);
fitStr.FDEF     = params(end-5);
fitStr.FDT      = params(end-4);
fitStr.FDW      = params(end-3);
fitStr.BGR      = params(end-2);
fitStr.BIN      = params(end-1);
fitStr.BCO      = params(end);
%% 4.4 - Storing the Eb(kx) dispersions
fitStr.nkx      = linspace(fitStr.kx_lims(1), fitStr.kx_lims(2), 1e3);
fitStr.neb      = {};
for i = 1:n
    fitStr.neb{i} = eb_calc(fitStr.nkx, fitStr.XLOC(i), fitStr.YLOC(i), fitStr.MSTAR(i));
end
end

%% DEFINING THE FUNCTION TO BE MINIMISED VIA GLOBAL OPTIMISATION METHODS
function MINFUN = minimize_function(x, ARPESObj)
    [ARPES, MODEL] = extract_arpes_and_model(x, ARPESObj.arpesStr, ARPESObj.fit_args.cTYPE, ARPESObj.fit_args.modelFunc);
    % - 1 - Extracting the DATA
    D   = ARPES.data;
    % - 2 - Extracting the MODEL
    M   = MODEL.data;
    % - 3 - Extracting the RESIDUALS (residuals = data - (model + background))
    R   = D - M;                  
    % - 5 - Extracting the CHI-SQUARED value to be minimised
    % MINFUN = std2(R.^2 ./ abs(M + B));      % This was found to yield the lowest value for R^2, rather than sum
    MINFUN = sum(sum(R.^2 ./ abs(M)));   % This is the chi-sqaured
end

%% DEFINING THE FUNCTION THAT EXTRACTS THE ARPES AND MODEL DATA CONSISTENTLY
function [ARPES, MODEL] = extract_arpes_and_model(x, arpesStr, cTYPE, modelFunc)
    % 1 - Initialising variables
    n = length(cTYPE);
    % - Extracting parabolic parameters
    BOFF    = x(1);                 % scalar of the band offset
    XLOC    = x(2:2+n-1);        	% scalar of the x-location of the min/max of the 2D parabolic ARPES dispersion.
    INT     = x(2+n:2+2*n-1);      	% scalar of the peak intensity of 2D PE curve.
    XFWHM   = x(2+2*n:2+3*n-1);    	% scalar of the x-axis FWHM for each Gaussian (k-resolution)
    YFWHM   = x(2+3*n:2+4*n-1);   	% scalar of the y-axis FWHM for each Gaussian (Eb-resolution)
    MSTAR   = x(2+4*n:2+5*n-1);   	% scalar of the effective mass, which governs the curvature of the parabola.
    % - Extracting the binding energy of each subband
    YLOC = [];
    for i = 1:n
        YLOC(i) = modelFunc{i}(BOFF);  	% scalar of the y-location of the min/max of the 2D parabolic ARPES dispersion.
    end
    % - Extracting background parameters
    FDEF    = x(end-5);         % scalar of the FDD Fermi-Level position.
    FDT     = x(end-4);         % scalar of the FDD temperature.
    FDW     = x(end-3);         % scalar of the FDD Gaussian width after convolution.
    BGR     = x(end-2);         % scalar of the gradient of the linear background.
    BIN     = x(end-1);         % scalar of the y-intercept of the linear background.
    BCO     = x(end);           % scalar of the constant background y-offset value.
    % - Extracting the kx and eb limits
    kx_lims = [min(arpesStr.kx(:)), max(arpesStr.kx(:))];
    eb_lims = [min(arpesStr.eb(:)), max(arpesStr.eb(:))];
    % 2 - Extracting the constructed MODEL data
    KX    = linspace(kx_lims(1), kx_lims(2), size(arpesStr.kx,2));
    EB    = linspace(eb_lims(1), eb_lims(2), size(arpesStr.eb,1))';
    MODEL_DATA    = ARPESCurve_FDDGpL(KX, EB, cTYPE, INT, XLOC, YLOC, XFWHM, YFWHM, MSTAR, FDEF, FDT, FDW, BGR, BIN, BCO);
%     MODEL_DATA    = ARPESCurveNorm_FDDGpL(KX, EB, cTYPE, INT, XLOC, YLOC, XFWHM, YFWHM, MSTAR, FDEF, FDT, FDW, BGR, BIN, BCO);
    ARPES_DATA    = arpesStr.data;
    % - Removing any NaN values
    MODEL_DATA(isnan(MODEL_DATA)) = 0;
    ARPES_DATA(isnan(ARPES_DATA)) = 0;
    % 3 - Defining the ARPES data structure
    ARPES               = struct();
    ARPES.kx           	= KX;
    ARPES.kx_step       = mean(diff(KX(:)));
    ARPES.kx_lims      	= kx_lims;
    ARPES.eb           	= EB;
    ARPES.eb_step       = mean(diff(EB(:)));
    ARPES.eb_lims     	= eb_lims;
    ARPES.data         	= ARPES_DATA;
    % 4 - Defining the MODEL data structure
    MODEL               = struct();
    MODEL.kx           	= KX;
    MODEL.kx_step      	= mean(diff(KX(:)));
    MODEL.kx_lims      	= kx_lims;
    MODEL.eb           	= EB;
    MODEL.eb_step     	= mean(diff(EB(:)));
    MODEL.eb_lims     	= eb_lims;
    MODEL.data        	= MODEL_DATA;
end