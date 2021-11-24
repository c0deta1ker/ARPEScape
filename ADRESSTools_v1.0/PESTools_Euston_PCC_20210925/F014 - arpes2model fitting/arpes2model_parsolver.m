function fitStr = arpes2model_parsolver(MODEL, ARPES, iparams, ibgrnd, solve_type)
% fitStr = arpes2model_parsolver(MODEL, ARPES, iparams, ibgrnd, solve_type)
%   Function that runs a a Global Optimisation Algorithm based on Simulated
%   Annealing to determine the best model fit to the input ARPES data and
%   model (initial conditions and bounds). The temperature of the Simulated
%   Annealing is set to walk a large span of the parameter space to ensure
%   the global minimum is obtained. The minimisation function is based on
%   minimising the standard deviation of the sum of the squared residuals.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   MODEL:          MATLAB data-structure that contains the model data.
%   -   ARPES:         	MATLAB data-structure that contains the ARPES data.
%   -   iparams:      	3 cells {x0}{lb}{ub} with 1x5 array: the model fit parameters [kxFWHM,ebFWHM,MFP,BOFF,INT]
%   -   ibgrnd:         3 cells {x0}{lb}{ub} with 1x3 vector of the background parameters: [MX,MY,C]
%   -   solve_type:     String of either "lsqnonlin", "simulannealbnd"
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information

%% Default parameters
% - Default input parameters to use
if nargin < 5; solve_type = "fmincon"; end
if nargin < 4
    ibgrnd{1} = [0, 0, 0];
    ibgrnd{2} = ibgrnd{1} - 0.10;
    ibgrnd{3} = ibgrnd{1} + 0.10;
end
if isempty(solve_type); solve_type = "fmincon"; end
if isempty(ibgrnd)
    ibgrnd{1} = [0, 0, 0];
    ibgrnd{2} = ibgrnd{1} - 0.10;
    ibgrnd{3} = ibgrnd{1} + 0.10;
end

%% - 1 - CONSTRAINING THE FIT PARAMETERS TO PHYSICAL VALUES
% -- Defining the initial conditions and limits of the parameters
x0      = [iparams{1}(:); ibgrnd{1}(:)];
lb      = [iparams{2}(:); ibgrnd{2}(:)];
ub      = [iparams{3}(:); ibgrnd{3}(:)];
% -- Applying constraints on the parameters, based on the data
% --- parameters
if x0(1) < 0; x0(1) = 0; end    % if kxFWHM < 0, then make it 0
if x0(2) < 0; x0(2) = 0; end    % if ebFWHM < 0, then make it 0
if x0(3) < min(MODEL.MFP(:)); x0(3) = min(MODEL.MFP(:)); end    % if MFP < min available, then make it the min
if x0(3) > max(MODEL.MFP(:)); x0(3) = max(MODEL.MFP(:)); end    % if MFP > max available, then make it the max
if x0(4) < min(MODEL.BOFF(:)); x0(4) = min(MODEL.BOFF(:)); end    % if BOFF < min available, then make it the min
if x0(4) > max(MODEL.BOFF(:)); x0(4) = max(MODEL.BOFF(:)); end    % if BOFF > max available, then make it the max
if x0(5) < 0; x0(5) = 0; end    % if INT < 0, then make it 0
% --- lower bounds
if lb(1) < 0; lb(1) = 0; end    % if kxFWHM < 0, then make it 0
if lb(2) < 0; lb(2) = 0; end    % if ebFWHM < 0, then make it 0
if lb(3) < min(MODEL.MFP(:)); lb(3) = min(MODEL.MFP(:)); end    % if MFP < min available, then make it the min
if lb(3) > max(MODEL.MFP(:)); lb(3) = max(MODEL.MFP(:)); end    % if MFP > max available, then make it the max
if lb(4) < min(MODEL.BOFF(:)); lb(4) = min(MODEL.BOFF(:)); end    % if BOFF < min available, then make it the min
if lb(4) > max(MODEL.BOFF(:)); lb(4) = max(MODEL.BOFF(:)); end    % if BOFF > max available, then make it the max
if lb(5) < 0; lb(5) = 0; end    % if INT < 0, then make it 0
% --- upper bounds
if ub(1) < 0; ub(1) = 0; end    % if kxFWHM < 0, then make it 0
if ub(2) < 0; ub(2) = 0; end    % if ebFWHM < 0, then make it 0
if ub(3) < min(MODEL.MFP(:)); ub(3) = min(MODEL.MFP(:)); end    % if MFP < min available, then make it the min
if ub(3) > max(MODEL.MFP(:)); ub(3) = max(MODEL.MFP(:)); end    % if MFP > max available, then make it the max
if ub(4) < min(MODEL.BOFF(:)); ub(4) = min(MODEL.BOFF(:)); end    % if BOFF < min available, then make it the min
if ub(4) > max(MODEL.BOFF(:)); ub(4) = max(MODEL.BOFF(:)); end    % if BOFF > max available, then make it the max
if ub(5) < 0; ub(5) = 0; end    % if INT < 0, then make it 0
% -- Determining the degrees of freedom for the fitting model
val	= lb + (ub-lb); DoF = 0;
for i = 1:length(val); if val(i) ~= lb(i); DoF = DoF + 1; end; end

%% - 2 - DEFINING THE ARPES DATA AND FITTING ARGUMENTS
% -- Defining global variable that stores all relevant model and data variables
ARPESObj        = struct();
ARPESObj.ARPES	= ARPES;
ARPESObj.MODEL	= MODEL;
% -- Appending the input arguments to the global variable
ARPESObj.fit_args.solve_type  = solve_type;
ARPESObj.fit_args.iparams   	= iparams;
ARPESObj.fit_args.ibgrnd    	= ibgrnd;

%% - 3 - RUNNING THE SOLVER
% -- Local solver properties
MaxFunEvals  = 1e5;
MaxIter     = 1e4;
FuncTol     = 1e-10;
StepTol     = 1e-10;
OptTol      = 1e-10;
ConTol      = 1e-10;
% -- Simulated annealing properties
TempFcn     = 'temperaturefast';
InitTemp    = 1e2;
ReannealInt = 1e2;
% -- Initialising parallel processing
nCPU_PCC_Work = 6;
nCPU_PSI_HPC_Merlin6 = 44;
myCluster = parcluster('local'); nCPU = myCluster.NumWorkers;
if isempty(gcp('nocreate')); poolobj = parpool('local', nCPU, 'IdleTimeout', Inf); end
%% - 3.1 - LOCAL SOLVER: FIND MINIMUM OF UNBOUNDED MULTIVARIATE FUNCTION
if solve_type == "fminunc"
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions('fminunc',...
        'MaxFunctionEvaluations', MaxFunEvals,...
        'FunctionTolerance', FuncTol,...
        'StepTolerance', StepTol,...
        'OptimalityTolerance', OptTol,...
        'MaxIterations', MaxIter,...
        'UseParallel', true);
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
        'ConstraintTolerance', ConTol,...
        'UseParallel', true);
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
        'ConstraintTolerance', ConTol,...
        'UseParallel', true);
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
% -- Shutting down parallel processing
% delete(poolobj);
%% - 4 - STORING THE FITS AS A MATLAB STRUCTURE
fitStr              = struct();
%% 4.1 - Storing the initial arguments
fitStr.solve_type	= solve_type;
fitStr.iparams      = iparams;
fitStr.ibgrnd       = ibgrnd;
%% 4.2 - Storing the ARPES data, background, model and minimisation variables
% -- Storing the domain and range
fitStr.kx           = MODEL.kx;
fitStr.eb           = MODEL.eb;
% -- x- and y-axis limits
fitStr.kx_lims      = [min(fitStr.kx(:)), max(fitStr.kx(:))];
fitStr.eb_lims      = [min(fitStr.eb(:)), max(fitStr.eb(:))];
% -- Extracting the  DATA
fitStr.D        = fit_data(params, ARPESObj);
% -- Extracting the MODEL
fitStr.M        = fit_model(params, ARPESObj);
% -- Extracting the BACKGROUND
fitStr.B        = fit_background(params, ARPESObj);
% -- Extracting the DATA - BACKGROUND
fitStr.DB       = fitStr.D - fitStr.B;
fitStr.MB     	= fitStr.M + fitStr.B;
% -- Extracting the RESIDUALS (R = D - (M + B))
fitStr.R        = fitStr.D - (fitStr.M + fitStr.B);
fitStr.Rx    	= mean(fitStr.R, 1);      % R in 1D along x
fitStr.Ry       = mean(fitStr.R, 2);      % R in 1D along y
% -- Extracting the CHI-SQUARED value to be minimised
fitStr.CHISQ2D	= fitStr.R.^2 ./ abs(fitStr.M + fitStr.B);
fitStr.CHISQx 	= mean(fitStr.CHISQ2D, 1);      % CHISQ in 1D along x
fitStr.CHISQy  	= mean(fitStr.CHISQ2D, 2);      % CHISQ in 1D along y
fitStr.CHISQ    = sum(sum(fitStr.CHISQ2D));
fitStr.MINFUN   = minimize_function(params, ARPESObj);
% -- Append the degrees of freedom
fitStr.DoF   	= DoF;
%% 4.3 - Storing the final fit variables for the ARPES model fit
fitStr.params	= params;
fitStr.kxFWHM  	= params(1);
fitStr.ebFWHM   = params(2);
fitStr.MFP      = params(3);
fitStr.BOFF     = params(4);
fitStr.INT      = params(5);
fitStr.MX       = params(6);
fitStr.MY       = params(7);
fitStr.C        = params(8);
end

%% DEFINING THE FUNCTION TO BE MINIMISED VIA GLOBAL OPTIMISATION METHODS
function MINFUN = minimize_function(x, ARPESObj)
    % - 1 - Extracting the DATA
    D   = fit_data(x, ARPESObj);
    % - 2 - Extracting the MODEL
    M   = fit_model(x, ARPESObj);
    % - 3 - Extracting the BACKGROUND
    B   = fit_background(x, ARPESObj);
    % - 4 - Extracting the RESIDUALS (residuals = data - (model + background))
    R   = D - (M + B);                  
    % - 5 - Extracting the CHI-SQUARED value to be minimised
    % MINFUN = std2(R.^2 ./ abs(M + B));      % This was found to yield the lowest value for R^2, rather than sum
    MINFUN = sum(sum(R.^2 ./ abs(M + B)));   % This is the chi-sqaured
end

%% DEFINING THE FUNCTION THAT DETERMINES THE ARPES DATA TO BE FITTED
function D = fit_data(x, ARPESObj)
    D   = ARPESObj.ARPES.data;
    D(isnan(D)) = 0;
end

%% DEFINING THE FUNCTION THAT DETERMINES THE MODEL
function M = fit_model(x, ARPESObj)
    % 1 - Initialising variables
    kxFWHM  = x(1);    % scalar of kx-resolution (Gaussian FWHM of k-broadening)
    ebFWHM  = x(2);    % scalar of eb-resolution (Gaussian FWHM of eb-broadening)
    MFP     = x(3);    % scalar of mean free path (associated with resolution of bands)
    BOFF    = x(4);    % scalar of band-offset
    INT     = x(5);    % scalar of the intensity scaling factor
    % 2 - Extracting the model curve based on the MFP and BOFF chosen
    [~, nMFP]       = min(abs(ARPESObj.MODEL.MFP - MFP));
    [~, nBOFF]      = min(abs(ARPESObj.MODEL.BOFF - BOFF));
    model_data   	= squeeze(ARPESObj.MODEL.data(nMFP,nBOFF,:,:))';
    % 3 - Gaussian broadening the model in x- and y-dimensions
    % -- Extracting the FWHM of the smoothing
    kxFWHM         = kxFWHM ./ abs(ARPESObj.MODEL.kx(2) - ARPESObj.MODEL.kx(1));
    ebFWHM         = ebFWHM ./ abs(ARPESObj.MODEL.eb(2) - ARPESObj.MODEL.eb(1));
    % -- Gaussian smoothing
    model_data   	= Gaco2(model_data, kxFWHM, ebFWHM);
    % 4 - Intensity scaling
    model_data      = INT .* model_data;
    % Assigning the final model
    M               = model_data;
    M(isnan(M))     = 0;
end

%% DEFINING THE FUNCTION THAT DETERMINES THE PLANAR BACKGROUND
function B = fit_background(x, ARPESObj)
    % 1 - Initialising variables
    MX  = x(end-2);   	% scalar of the x-gradient of the plane.
    MY  = x(end-1);    	% scalar of the y-gradient of the plane.
    C   = x(end);     	% scalar of the constant off-set of the plane.
    % 2 - Creating a planar background to be subtracted
    B  	= MX * ARPESObj.MODEL.kx + MY * ARPESObj.MODEL.eb + C;
    B(isnan(B)) = 0;
end