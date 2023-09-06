function fitStr = pes2boff_solver(pesStr, modelStr, cTYPE, iparams, bTYPE, ibgrnd, solve_type)
% fitStr = pes2boff_solver(pesStr, modelStr, cTYPE, iparams, bTYPE, ibgrnd, solve_type)
%   Function that runs an optimisation algorithm to solve for the best fit
%   to the PES data using N curves that are integrated over a predefined
%   potential profile that is defined by a particular model data structure and
%   modulated by the mean free path.
%   
%   REQ. FUNCTIONS: none
%   
%   IN:
%   -   pesStr:         data structure that contains the PES data.
%   -   modelStr:       data structure that contains the model data.
%   -   cTYPE:          1xN vector of the type of curve to use for the nth state.
%   -   iparams:   	    3 cells {x0}{lb}{ub} with Nx8 arrays: the n'th peak parameters [BE,INT,FWHM,MR,LSE,LSI,LSW,ASY]
%   -   bTYPE:          string of the type of background to use for fitting. Default: "Poly" ("none", "Poly", "Shir", "LinShir", "StepFDDGpL", "StepFDDGsL")
%   -   ibgrnd:         4 cells {x0}{lb}{ub}{args} with 1x3 vectors of the background parameters: x0=lb=ub=[LHS,RHS,BGR], or argument of the background type args = []
%   -   solve_type:     string of the optimisation algorithm to use; "lsqcurvefit" (very fast, good), "lsqnonlin" (fast, not good), "simulannealbnd" (slow, best). 
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information

%% Default parameters
% - Default input parameters to use
if nargin < 6; solve_type = "lsqcurvefit"; end
if nargin < 5
    ibgrnd{1} = [...
        mean(pesStr.xdat(:)) - abs(0.25*range(pesStr.xdat(:))),...
        mean(pesStr.xdat(:)) + abs(0.25*range(pesStr.xdat(:))),...
        0];
    ibgrnd{2} = [0, 0, -abs(0.05*range(pesStr.ydat(:)))];
    ibgrnd{3} = [0, 0, +abs(0.05*range(pesStr.ydat(:)))];
    ibgrnd{4} = {1};
end
if nargin < 4; bTYPE = "Poly"; end
if isempty(solve_type); solve_type = "lsqcurvefit"; end
if isempty(ibgrnd)
    ibgrnd{1} = [...
        mean(pesStr.xdat(:)) - abs(0.25*range(pesStr.xdat(:))),...
        mean(pesStr.xdat(:)) + abs(0.25*range(pesStr.xdat(:))),...
        0];
    ibgrnd{2} = [0, 0, -abs(0.05*range(pesStr.ydat(:)))];
    ibgrnd{3} = [0, 0, +abs(0.05*range(pesStr.ydat(:)))];
    ibgrnd{4} = {1};
end
if isempty(bTYPE); bTYPE = "Poly"; end

%% - 1 - CONSTRAINING THE FIT PARAMETERS TO PHYSICAL VALUES
% -- Total number of states to be fitted
nSTATES	= length(cTYPE);
% -- Defining the initial conditions and limits of the optimisation method
x0 = [iparams{1}(:); ibgrnd{1}(:)];
lb = [iparams{2}(:); ibgrnd{2}(:)];
ub = [iparams{3}(:); ibgrnd{3}(:)];
% -- Defining the initial conditions and limits of the optimisation method
for i = 1:nSTATES
    indx = (2+i):nSTATES:(length(x0)-3); % +2 to account for MFP and BOFF                    
    X0 = x0(indx); LB = lb(indx); UB = ub(indx);
    % --- parameters
    if X0(2) < 0; x0(indx(2)) = 0; end    % if INT < 0, then make it 0
    if X0(3) < 0; x0(indx(3)) = 0; end    % if FWHM < 0, then make it 0
    if X0(4) < 0; x0(indx(4)) = 0; end    % if MR < 0, then make it 0
    if X0(4) > 1; x0(indx(4)) = 1; end    % if MR > 1, then make it 1
    if X0(6) < 0; x0(indx(6)) = 0; end    % if LSI < 0, then make it 0
    if X0(7) < 0; x0(indx(7)) = 0; end    % if LSW < 0, then make it 0
    if X0(8) < 0; x0(indx(8)) = 0; end    % if ASY < 0, then make it 0
    % --- lower bounds
    if LB(2) < 0; lb(indx(2)) = 0; end    % if INT < 0, then make it 0
    if LB(3) < 0; lb(indx(3)) = 0; end    % if FWHM < 0, then make it 0
    if LB(4) < 0; lb(indx(4)) = 0; end    % if MR < 0, then make it 0
    if LB(4) > 1; lb(indx(4)) = 1; end    % if MR > 1, then make it 1
    if LB(6) < 0; lb(indx(6)) = 0; end    % if LSI < 0, then make it 0
    if LB(7) < 0; lb(indx(7)) = 0; end    % if LSW < 0, then make it 0
    if LB(8) < 0; lb(indx(8)) = 0; end    % if ASY < 0, then make it 0
    % --- lower bounds
    if UB(2) < 0; ub(indx(2)) = 0; end    % if INT < 0, then make it 0
    if UB(3) < 0; ub(indx(3)) = 0; end    % if FWHM < 0, then make it 0
    if UB(4) < 0; ub(indx(4)) = 0; end    % if MR < 0, then make it 0
    if UB(4) > 1; ub(indx(4)) = 1; end    % if MR > 1, then make it 1
    if UB(6) < 0; ub(indx(6)) = 0; end    % if LSI < 0, then make it 0
    if UB(7) < 0; ub(indx(7)) = 0; end    % if LSW < 0, then make it 0
    if UB(8) < 0; ub(indx(8)) = 0; end    % if ASY < 0, then make it 0
    % --- ensuring the fit window is within the domain
    if x0(end-2) < min(pesStr.xdat(:)); x0(end-2) = min(pesStr.xdat(:)); end
    if x0(end-1) > max(pesStr.xdat(:)); x0(end-1) = max(pesStr.xdat(:)); end 
end
% -- Determining the degrees of freedom for the fitting model
val	= lb + (ub-lb); DoF = 0;
for i = 1:length(val); if val(i) ~= lb(i); DoF = DoF + 1; end; end

%% - 2 - DEFINING THE XPS DATA AND FITTING ARGUMENTS
% -- Extracting the data to be fitted
[roi_xdat, roi_ydat, roi_bgrnd] = PESBackground(pesStr.xdat, pesStr.ydat, bTYPE, x0(end-2), x0(end-1), x0(end), ibgrnd{4});
% -- Defining a structure that stores all relevant model and data variables
XPSObj              = pesStr;
XPSObj.roi_xdat     = roi_xdat;
XPSObj.roi_ydat     = roi_ydat;
XPSObj.roi_bgrnd	= roi_bgrnd;
XPSObj.roi_ydat_sbtr	= XPSObj.roi_ydat - XPSObj.roi_bgrnd;
% -- Appending the input arguments to the global variable
XPSObj.fit_args.solve_type  = solve_type;
XPSObj.fit_args.nSTATES   	= length(cTYPE);
XPSObj.fit_args.cTYPE    	= cTYPE;
XPSObj.fit_args.iparams     = iparams;
XPSObj.fit_args.bTYPE    	= bTYPE;
XPSObj.fit_args.ibgrnd      = ibgrnd;
XPSObj.modelStr             = modelStr;

%% - 3 - RUNNING THE SOLVER
% -- Local solver properties
MaxFunEvals = 1e5;
MaxIter     = 5e3;
FuncTol     = 1e-7; % 1e-7;
StepTol     = 1e-7; % 1e-7;
OptTol      = 1e-7; % 1e-7;
ConTol      = 1e-7; % 1e-7;
FinDiffRelStep = 1e-5; % 1e-5;
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
    [params,rnorm,resid,exitflag,output,lambda,jacobian]  = lsqcurvefit(@(x,roi_xdat) full_pes_function(x,roi_xdat,XPSObj), x0, roi_xdat, roi_ydat, lb, ub, options);  
%% - 3.2 - LOCAL SOLVER: UNBOUNDED LEAST SQUARES NONLINEAR REGRESSION FIT
elseif solve_type == "nlinfit"
    % -- Defining the optimisation options for the simulated annealing method
    options = statset(...
        'MaxFunEvals', MaxFunEvals,...
        'MaxIter', MaxIter,...
        'TolFun', FuncTol);
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(roi_xdat, roi_ydat, @(x,roi_xdat) full_pes_function(x,roi_xdat,XPSObj), x0, options);
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
    [params,fval,exitflag,output,grad,hessian] = fminunc(@(x) minimize_function(x,XPSObj), x0, options); 
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
    [params,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x) minimize_function(x,XPSObj), x0, [], [], [], [], lb, ub, [], options);  
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
    [params,rnorm,resid,exitflag,output,lambda,jacobian] = lsqnonlin(@(x) minimize_function(x,XPSObj), x0, lb, ub, options);    
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
    [params,fval,exitflag,output] = simulannealbnd(@(x) minimize_function(x,XPSObj), x0, lb, ub, options);
end

%% - 4 - STORING THE FITS AS A MATLAB STRUCTURE
fitStr              = struct();
%% 4.1 - Storing the initial arguments
fitStr.solve_type	= solve_type;
fitStr.cTYPE        = cTYPE;
fitStr.iparams      = iparams;
fitStr.bTYPE        = bTYPE;
fitStr.ibgrnd       = ibgrnd;
fitStr.nSTATES      = length(cTYPE);
%% 4.2 - Storing the XPS data, background, model and minimisation variables
% -- Extracting the original data 
fitStr.xdat     = pesStr.xdat;
fitStr.ydat 	= pesStr.ydat;
% -- Extracting the X-DOMAIN and DATA
[fitStr.X, fitStr.D] = fit_data(params, XPSObj);
% -- Extracting the MODEL
fitStr.M        = fit_model(params, XPSObj);
% -- Extracting the BACKGROUND
fitStr.B        = fit_background(params, XPSObj);
% -- Extracting the DATA - BACKGROUND
fitStr.DB       = fitStr.D - fitStr.B;
fitStr.MB       = fitStr.M + fitStr.B;
% -- Extracting the RESIDUALS (R = D - (M + B))
fitStr.R        = fitStr.D - (fitStr.M + fitStr.B);
% -- Extracting the CHI-SQUARED value to be minimised
fitStr.CHISQ    = sum(fitStr.R.^2 ./ abs(fitStr.M + fitStr.B));
fitStr.MINFUN   = minimize_function(params, XPSObj);
% -- Append the degrees of freedom
fitStr.DoF   	= DoF;
%% 4.3 - Storing the final fit variables for the PES background
fitStr.bPARAMS	= params(end-2:end);
fitStr.LHS      = fitStr.bPARAMS(1);
fitStr.RHS      = fitStr.bPARAMS(2);
fitStr.BGR      = fitStr.bPARAMS(3);
%% 4.4 - Storing the final fit variables of each PES curve component
fitStr.MFP  	= params(1);
fitStr.BOFF 	= params(2);
fitStr.XX     	= linspace(min(fitStr.X(:)), max(fitStr.X(:)), 5e3)';
fitStr.YY    	= zeros(size(fitStr.XX));
% -- Storing the potential profile
imodelStr      	= mstheory_interp_spectra(modelStr, fitStr.BOFF, fitStr.MFP);
fitStr.ZPOT 	= linspace(min(imodelStr.ZPOT(:)), max(imodelStr.ZPOT(:)), 1e3);
fitStr.EPOT     = interp1(imodelStr.ZPOT, imodelStr.EPOT, fitStr.ZPOT);
for i = 1:fitStr.nSTATES
    % --- Best fit components on new domain
    fitStr.cPARAMS(i,:) = params((2+i):fitStr.nSTATES:(end-3));
   	[fitStr.cYY(:,i), fitStr.potYY(:,:,i)] = PESCurve_POT(fitStr.XX, fitStr.cTYPE(i),...
        fitStr.cPARAMS(i,1), fitStr.cPARAMS(i,2), fitStr.cPARAMS(i,3),...
        fitStr.cPARAMS(i,4), fitStr.cPARAMS(i,5), fitStr.cPARAMS(i,6),...
        fitStr.cPARAMS(i,7), fitStr.cPARAMS(i,8),...
        fitStr.MFP, fitStr.ZPOT, fitStr.EPOT);
    fitStr.YY = fitStr.YY + fitStr.cYY(:,i);
    % --- Storing each curve parameter
    fitStr.DCL(i)   = fitStr.cPARAMS(i,1);
    fitStr.INT(i)   = fitStr.cPARAMS(i,2);
    fitStr.FWHM(i)  = fitStr.cPARAMS(i,3);
    fitStr.MR(i)    = fitStr.cPARAMS(i,4);
    fitStr.LSE(i)   = fitStr.cPARAMS(i,5);
    fitStr.LSI(i)   = fitStr.cPARAMS(i,6);
    fitStr.LSW(i)   = fitStr.cPARAMS(i,7);
    fitStr.ASY(i)   = fitStr.cPARAMS(i,8);
    % --- Storing energy information 
    fitStr.BE_Z0(i)	= fitStr.DCL(i) - fitStr.BOFF;
    [~,I] = max(fitStr.cYY(:,i)); fitStr.BE_MAX(i)= fitStr.XX(I);
    % --- Storing the curve area information
    fitStr.AREA(i)	= trapz(fitStr.XX, fitStr.cYY(:,i));
end
% --- Storing the normalised area contribution for quantification
fitStr.AREA0        = fitStr.AREA ./ trapz(fitStr.XX, fitStr.YY);
% - Extracting the total potential profile of all components
fitStr.potYY_tot         = sum(fitStr.potYY, 3);
fitStr.potYY_comp        = squeeze(sum(fitStr.potYY, 2));
for i = 1:nSTATES
    fitStr.potYY_comp(:,i) = fitStr.INT(i) .* (fitStr.potYY_comp(:,i) ./ max(fitStr.potYY_comp(:,i)));
end

end

%% DEFINING THE RESIDUAL MINIMISATION FUNCTION TO BE MINIMISED
function MINFUN = minimize_function(x, XPSObj)
    % - 1 - Extracting the DATA
    [~, D] = fit_data(x, XPSObj);
    % - 2 - Extracting the MODEL
    M = fit_model(x, XPSObj);
    % - 3 - Extracting the BACKGROUND
    B = fit_background(x, XPSObj);
    % - 4 - Extracting the RESIDUALS (residuals = data - (model + background))
    R   = D - (M + B);
    % - 5 - Extracting the CHI-SQUARED value to be minimised
    MINFUN = sum(R.^2 ./ abs(M + B));
%     MINFUN = std(R.^2 ./ abs(M + B));
end

%% DEFINING THE FUNCTION THAT EXTRACTS THE XPS DATA TO BE FITTED
function [X, D] = fit_data(x, XPSObj)
    [X, D, ~] = PESBackground(XPSObj.xdat, XPSObj.ydat, XPSObj.fit_args.bTYPE, x(end-2), x(end-1), x(end), XPSObj.fit_args.ibgrnd{4});
    D(isnan(D)) = 0;
    if size(D, 2) > 1; D = D'; end
end

%% DEFINING THE FUNCTION THAT DETERMINES THE TOTAL PES BACKGROUND
function B = fit_background(x, XPSObj)
    [~, ~, B] = PESBackground(XPSObj.xdat, XPSObj.ydat, XPSObj.fit_args.bTYPE, x(end-2), x(end-1), x(end), XPSObj.fit_args.ibgrnd{4});
    B(isnan(B)) = 0;
    if size(B, 2) > 1; B = B'; end
end

%% DEFINING THE FUNCTION THAT DETERMINES THE TOTAL PES CURVE FIT
function M = fit_model(x, XPSObj)
    % - 1 - Extracting the potential profile
    imodelStr     	= mstheory_interp_spectra(XPSObj.modelStr, x(2), x(1));
    ZPOT            = linspace(min(imodelStr.ZPOT(:)), max(imodelStr.ZPOT(:)), 1e3);
    EPOT            = interp1(imodelStr.ZPOT, imodelStr.EPOT, ZPOT);
    % - 2 - Extract the curve component for each state
    comp_int = {};
    for i = 1:XPSObj.fit_args.nSTATES
        % -- Extracting the arguments for the component curve
        pes_args    = x((2+i):XPSObj.fit_args.nSTATES:(end-3));
        % -- Extracting the component intensities
        comp_int{i} = PESCurve_POT(XPSObj.roi_xdat, XPSObj.fit_args.cTYPE(i),...
            pes_args(1), pes_args(2), pes_args(3),...
            pes_args(4), pes_args(5), pes_args(6),...
            pes_args(7), pes_args(8),...
            x(1), ZPOT, EPOT);
    end
    % - 3 - Sum up all of the curve components to get the final model PES curve
    pes_int = zeros(size(XPSObj.roi_xdat));
    for i = 1:XPSObj.fit_args.nSTATES
        pes_int = pes_int + comp_int{i};
    end
    M = pes_int;
    M(isnan(M)) = 0;
    if size(M, 2) > 1; M = M'; end
end

%% DEFINING THE FINAL FUNCTION THAT DESCRIBES THE PES CURVE TO BE FITTED
function MB = full_pes_function(x, xdat, XPSObj)
    % - 1 - Extracting the potential profile
    imodelStr     	= mstheory_interp_spectra(XPSObj.modelStr, x(2), x(1));
    ZPOT            = linspace(min(imodelStr.ZPOT(:)), max(imodelStr.ZPOT(:)), 1e3);
    EPOT            = interp1(imodelStr.ZPOT, imodelStr.EPOT, ZPOT);
    % - 2 - Extract the curve component for each state
    comp_int = {};
    for i = 1:XPSObj.fit_args.nSTATES
        % -- Extracting the arguments for the component curve
        pes_args    = x((2+i):XPSObj.fit_args.nSTATES:(end-3));
        % -- Extracting the component intensities
        comp_int{i} = PESCurve_POT(xdat, XPSObj.fit_args.cTYPE(i),...
            pes_args(1), pes_args(2), pes_args(3),...
            pes_args(4), pes_args(5), pes_args(6),...
            pes_args(7), pes_args(8),...
            x(1), ZPOT, EPOT);
    end
    % - 3 - Sum up all of the curve components to get the final model PES curve
    pes_int = zeros(size(xdat));
    for i = 1:XPSObj.fit_args.nSTATES
        pes_int = pes_int + comp_int{i};
    end
    M = pes_int;
    % - 4 - Determine the background to be used
    [B_xdat, ~, B] = PESBackground(xdat, XPSObj.roi_ydat, XPSObj.fit_args.bTYPE, x(end-2), x(end-1), x(end), XPSObj.fit_args.ibgrnd{4});
    B = interp1(B_xdat, B, xdat);
    % - 5 - Final output is the sum of the model and background
    MB = M + B;
    MB(isnan(MB)) = 0; if size(MB, 2) > 1; MB = MB'; end
end