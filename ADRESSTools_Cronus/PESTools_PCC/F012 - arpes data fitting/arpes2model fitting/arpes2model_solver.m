function fitStr = arpes2model_solver(arpesStr, modelStr, iparams, ibgrnd, solve_type, kx_lims, eb_lims)
% fitStr = arpes2model_solver(arpesStr, modelStr, iparams, ibgrnd, solve_type, kx_lims, eb_lims)
%   Function that runs a a Global Optimisation Algorithm based on Simulated
%   Annealing to determine the best model fit to the input ARPES data and
%   model (initial conditions and bounds). The temperature of the Simulated
%   Annealing is set to walk a large span of the parameter space to ensure
%   the global minimum is obtained. The minimisation function is based on
%   minimising the value of chi-squared.
%
%   IN:
%   -   arpesStr:       MATLAB data-structure that contains the initial ARPES data.
%   -   modelStr:       MATLAB data-structure that contains the initial model data.
%   -   iparams:      	3 cells {x0}{lb}{ub} with [1×5] array: the model fit parameters [K0,KFWHM,EFWHM,MFP,BOFF,INT]
%   -   ibgrnd:         3 cells {x0}{lb}{ub} with [1×3] vector of the background parameters: [MX,MY,C]
%   -   solve_type:     String that sets the optimisation routine to use:
%                           "fminunc"
%                           "fmincon"
%                           "lsqnonlin"
%                           "simulannealbnd"
%   -   kx_lims:      	[1×2] vector of [minKx, maxKx], which defines the consistent kx fit window
%   -   eb_lims:        [1×2] vector of [minEb, maxEb], which defines the consistent eb fit window
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information

%% Default parameters & validity checks
% - Default input parameters to use
if nargin < 7; eb_lims = mean(modelStr.eb(:)) + [-1,1].*0.95*range(modelStr.eb(:)); end
if nargin < 6; kx_lims = mean(modelStr.kx(:)) + [-1,1].*0.95*range(modelStr.kx(:)); end
if nargin < 5; solve_type = "fmincon"; end
if nargin < 4; ibgrnd{1} = [0, 0, 0]; ibgrnd{2} = ibgrnd{1} - 0.01; ibgrnd{3} = ibgrnd{1} + 0.01; end
if isempty(eb_lims); eb_lims = mean(modelStr.eb(:)) + [-1,1].*0.95*range(modelStr.eb(:)); end
if isempty(kx_lims); kx_lims = mean(modelStr.eb(:)) + [-1,1].*0.95*range(modelStr.kx(:)); end
if isempty(solve_type); solve_type = "fmincon"; end
if isempty(ibgrnd); ibgrnd{1} = [0, 0, 0]; ibgrnd{2} = ibgrnd{1}-0.01; ibgrnd{3} = ibgrnd{1}+0.01; end
% -- Validity check on the input variables
if kx_lims(1) < 0.95*min(modelStr.kx(:)); kx_lims(1) = 0.95*min(modelStr.kx(:)); end
if kx_lims(2) > 0.95*max(modelStr.kx(:)); kx_lims(2) = 0.95*max(modelStr.kx(:)); end
if eb_lims(1) < 0.95*min(modelStr.eb(:)); eb_lims(1) = 0.95*min(modelStr.eb(:)); end
if eb_lims(2) > 0.95*max(modelStr.eb(:)); eb_lims(2) = 0.95*max(modelStr.eb(:)); end

%% - 1 - CONSTRAINING THE FIT PARAMETERS TO PHYSICAL VALUES
% -- Defining the initial conditions and limits of the parameters
x0      = [iparams{1}(:); ibgrnd{1}(:)];
lb      = [iparams{2}(:); ibgrnd{2}(:)];
ub      = [iparams{3}(:); ibgrnd{3}(:)];
% -- Applying constraints on the parameters, based on the data
% --- parameters
if x0(2) < 0; x0(2) = 0; end                                            % if kxFWHM < 0, then make it 0
if x0(3) < 0; x0(3) = 0; end                                            % if ebFWHM < 0, then make it 0
if x0(4) < min(modelStr.MFP(:)); x0(4) = min(modelStr.MFP(:)); end      % if MFP < min available, then make it the min
if x0(4) > max(modelStr.MFP(:)); x0(4) = max(modelStr.MFP(:)); end      % if MFP > max available, then make it the max
if x0(5) < min(modelStr.BOFF(:)); x0(5) = min(modelStr.BOFF(:)); end  	% if BOFF < min available, then make it the min
if x0(5) > max(modelStr.BOFF(:)); x0(5) = max(modelStr.BOFF(:)); end  	% if BOFF > max available, then make it the max
if x0(6) < 0; x0(6) = 0; end                                            % if INT < 0, then make it 0
% --- lower bounds
if lb(2) < 0; lb(2) = 0; end                                            % if kxFWHM < 0, then make it 0
if lb(3) < 0; lb(3) = 0; end                                            % if ebFWHM < 0, then make it 0
if lb(4) < min(modelStr.MFP(:)); lb(4) = min(modelStr.MFP(:)); end   	% if MFP < min available, then make it the min
if lb(4) > max(modelStr.MFP(:)); lb(4) = max(modelStr.MFP(:)); end    	% if MFP > max available, then make it the max
if lb(5) < min(modelStr.BOFF(:)); lb(5) = min(modelStr.BOFF(:)); end   	% if BOFF < min available, then make it the min
if lb(5) > max(modelStr.BOFF(:)); lb(5) = max(modelStr.BOFF(:)); end  	% if BOFF > max available, then make it the max
if lb(6) < 0; lb(6) = 0; end                                            % if INT < 0, then make it 0
% --- upper bounds
if ub(2) < 0; ub(2) = 0; end                                            % if kxFWHM < 0, then make it 0
if ub(3) < 0; ub(3) = 0; end                                            % if ebFWHM < 0, then make it 0
if ub(4) < min(modelStr.MFP(:)); ub(4) = min(modelStr.MFP(:)); end      % if MFP < min available, then make it the min
if ub(4) > max(modelStr.MFP(:)); ub(4) = max(modelStr.MFP(:)); end      % if MFP > max available, then make it the max
if ub(5) < min(modelStr.BOFF(:)); ub(5) = min(modelStr.BOFF(:)); end    % if BOFF < min available, then make it the min
if ub(5) > max(modelStr.BOFF(:)); ub(5) = max(modelStr.BOFF(:)); end    % if BOFF > max available, then make it the max
if ub(6) < 0; ub(6) = 0; end                                            % if INT < 0, then make it 0
% -- Determining the degrees of freedom for the fitting model
val	= lb + (ub-lb); DoF = 0;
for i = 1:length(val); if val(i) ~= lb(i); DoF = DoF + 1; end; end

%% - 2 - DEFINING THE ARPES DATA AND FITTING ARGUMENTS
% -- Defining global variable that stores all relevant model and data variables
ARPESObj                        = struct();
ARPESObj.arpesStr               = arpesStr;
ARPESObj.modelStr               = modelStr;
% -- Appending the input arguments to the global variable
ARPESObj.fit_args.solve_type    = solve_type;
ARPESObj.fit_args.iparams   	= iparams;
ARPESObj.fit_args.ibgrnd    	= ibgrnd;
ARPESObj.fit_args.kx_lims   	= kx_lims;
ARPESObj.fit_args.eb_lims    	= eb_lims;

%% - 3 - RUNNING THE SOLVER
% -- Local solver properties
MaxFunEvals = 1e5;
MaxIter     = 1e4;
FuncTol     = 1e-12;
StepTol     = 1e-12;
OptTol      = 1e-12;
ConTol      = 1e-12;
FinDiffRelStep = 1e-5;
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
    [params,fval,exitflag,output,grad,hessian] = fminunc(@(x) minimize_function(x,ARPESObj), x0, options); 
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
        'FinDiffRelStep', FinDiffRelStep,...
        'FinDiffType', 'central');
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
fitStr.iparams      = iparams;
fitStr.ibgrnd       = ibgrnd;
%% 4.2 - Storing the ARPES data, background, model and minimisation variables
[ARPES, MODEL]      = extract_arpes_and_model(params, ARPESObj.arpesStr, ARPESObj.modelStr, kx_lims, eb_lims);
% -- Storing the domain and range
fitStr.kx           = MODEL.kx;
fitStr.eb           = MODEL.eb;
% -- x- and y-axis limits
fitStr.kx_lims      = [min(fitStr.kx(:)), max(fitStr.kx(:))];
fitStr.eb_lims      = [min(fitStr.eb(:)), max(fitStr.eb(:))];
% -- Extracting the  DATA
fitStr.D        = fit_data(params, ARPES);
% -- Extracting the MODEL
fitStr.M        = fit_model(params, MODEL);
% -- Extracting the BACKGROUND
fitStr.B        = fit_background(params, MODEL);
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
fitStr.K0       = params(1);
fitStr.KFWHM  	= params(2);
fitStr.EFWHM    = params(3);
fitStr.MFP      = params(4);
fitStr.BOFF     = params(5);
fitStr.INT      = params(6);
fitStr.MX       = params(7);
fitStr.MY       = params(8);
fitStr.C        = params(9);
end

%% DEFINING THE FUNCTION TO BE MINIMISED VIA GLOBAL OPTIMISATION METHODS
function MINFUN = minimize_function(x, ARPESObj)
    [ARPES, MODEL] = extract_arpes_and_model(x, ARPESObj.arpesStr, ARPESObj.modelStr, ARPESObj.fit_args.kx_lims, ARPESObj.fit_args.eb_lims);
    % - 1 - Extracting the DATA
    D   = fit_data(x, ARPES);
    % - 2 - Extracting the MODEL
    M   = fit_model(x, MODEL);
    % - 3 - Extracting the BACKGROUND
    B   = fit_background(x, MODEL);
    % - 4 - Extracting the RESIDUALS (residuals = data - (model + background))
    R   = D - (M + B);                  
    % - 5 - Extracting the CHI-SQUARED value to be minimised
    MINFUN = sum(sum(R.^2 ./ abs(M + B)));   % This is the chi-sqaured
end

%% DEFINING THE FUNCTION TO BE MINIMISED VIA GLOBAL OPTIMISATION METHODS
function MINFUN = minimize_function_std2(x, ARPESObj)
    [ARPES, MODEL] = extract_arpes_and_model(x, ARPESObj.arpesStr, ARPESObj.modelStr, ARPESObj.fit_args.kx_lims, ARPESObj.fit_args.eb_lims);
    % - 1 - Extracting the DATA
    D   = fit_data(x, ARPES); 
    % - 2 - Extracting the MODEL
    M   = fit_model(x, MODEL); 
    % - 3 - Extracting the BACKGROUND
    B   = fit_background(x, MODEL);
    % - 4 - Extracting the RESIDUALS (residuals = data - (model + background))
    R   = D - (M + B);                  
    % - 5 - Extracting the CHI-SQUARED value to be minimised
    MINFUN = std2(R.^2 ./ abs(M + B));      % This was found to yield the lowest value for R^2, rather than sum
end

%% DEFINING THE FUNCTION THAT DETERMINES THE ARPES DATA TO BE FITTED
function D = fit_data(x, ARPES)
    D   = ARPES.data;
    D(isnan(D)) = 0;
end

%% DEFINING THE FUNCTION THAT DETERMINES THE MODEL
function M = fit_model(x, MODEL)
    % 1 - Initialising variables
    kxFWHM  = x(2);    % scalar of kx-resolution (Gaussian FWHM of k-broadening)
    ebFWHM  = x(3);    % scalar of eb-resolution (Gaussian FWHM of eb-broadening)
    INT     = x(6);    % scalar of the intensity scaling factor
    % -- Constraints
    if kxFWHM < 0; kxFWHM = 0; end
    if ebFWHM < 0; ebFWHM = 0; end
    if INT < 0; INT = 0; end
    % 2 - Extracting the model curve based on the MFP and BOFF chosen
    model_data   	= MODEL.data;
    % 3 - Gaussian broadening the model in x- and y-dimensions
    % -- Extracting the FWHM of the smoothing
    kxFWHM         = kxFWHM ./ abs(MODEL.kx_step);
    ebFWHM         = ebFWHM ./ abs(MODEL.eb_step);
    % -- Gaussian smoothing
    model_data   	= Gaco2(model_data, kxFWHM, ebFWHM);
    model_data      = INT .* model_data;
    % Assigning the final model
    M               = model_data;
    M(isnan(M))     = 0;
end

%% DEFINING THE FUNCTION THAT DETERMINES THE PLANAR BACKGROUND
function B = fit_background(x, MODEL)
    % 1 - Initialising variables
    MX  = x(end-2);   	% scalar of the x-gradient of the plane.
    MY  = x(end-1);    	% scalar of the y-gradient of the plane.
    C   = x(end);     	% scalar of the constant off-set of the plane.
    % 2 - Creating a planar background to be subtracted
    B  	= MX * MODEL.kx + MY * MODEL.eb + C;
    B(isnan(B)) = 0;
end

%% DEFINING THE FUNCTION THAT EXTRACTS THE ARPES AND MODEL DATA CONSISTENTLY
function [ARPES, MODEL] = extract_arpes_and_model(x, arpesStr, modelStr, kx_lims, eb_lims)
    % 1 - Initialising variables
    K0      = x(1);    % scalar of kx-origin offset
    MFP     = x(4);    % scalar of mean free path (associated with resolution of bands)
    BOFF    = x(5);    % scalar of band-offset
    % 2 - Cropping the ARPES data to extract dimensional size
    [xDat_crop, yDat_crop, ~] = data_crop2D(arpesStr.kx, arpesStr.eb, arpesStr.data, kx_lims, eb_lims);
    kx_length   = size(xDat_crop,2);
    eb_length   = size(yDat_crop,1);
    % 3 - Extracting the interpolated MODEL data
    imodelStr	= mstheory_interp_spectra(modelStr, BOFF, MFP);
    % 4 - Making the ARPES and MODEL data have a consistent domain
    KX          = linspace(kx_lims(1), kx_lims(2), kx_length);
    EB          = linspace(eb_lims(1), eb_lims(2), eb_length)';
    ARPES_DATA  = interpn(arpesStr.eb(:,1), arpesStr.kx(1,:)-K0, arpesStr.data, EB, KX);
    MODEL_DATA  = interpn(imodelStr.eb, imodelStr.kx, imodelStr.data, EB, KX);
    % 5 - Data renormalisation
    ARPES_DATA = ARPES_DATA - min(ARPES_DATA(:)); ARPES_DATA = ARPES_DATA ./ max(ARPES_DATA(:));
    MODEL_DATA = MODEL_DATA - min(MODEL_DATA(:)); MODEL_DATA = MODEL_DATA ./ max(MODEL_DATA(:));
    ARPES_DATA(isnan(ARPES_DATA)) = 0;
    MODEL_DATA(isnan(MODEL_DATA)) = 0;
    % 6 - Defining the ARPES data structure
    ARPES.kx           	= KX;
    ARPES.kx_step       = mean(diff(KX(:)));
    ARPES.kx_lims      	= kx_lims;
    ARPES.eb           	= EB;
    ARPES.eb_step       = mean(diff(EB(:)));
    ARPES.eb_lims     	= eb_lims;
    ARPES.data         	= ARPES_DATA;
    % 7 - Defining the MODEL data structure
    MODEL.MFP          	= MFP;
    MODEL.BOFF        	= BOFF;
    MODEL.kx           	= KX;
    MODEL.kx_step      	= mean(diff(KX(:)));
    MODEL.kx_lims      	= kx_lims;
    MODEL.eb           	= EB;
    MODEL.eb_step     	= mean(diff(EB(:)));
    MODEL.eb_lims     	= eb_lims;
    MODEL.data        	= MODEL_DATA;
end