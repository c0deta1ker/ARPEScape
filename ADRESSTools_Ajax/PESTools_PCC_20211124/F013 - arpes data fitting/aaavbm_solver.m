function fitStr = xps_solver(xpsStr, cTYPE, iPESCurves, dPESCurves, bTYPE, iPESBgrnd, dPESBgrnd, solve_type)
% fitStr = xps_solver(xpsStr, cTYPE, iPESCurves, dPESCurves, bTYPE, iPESBgrnd, dPESBgrnd, solve_type)
%   Function that runs a a Global Optimisation Algorithm based on Simulated
%   Annealing to determine the best curve fit to the input XPS data and
%   model (initial conditions and bounds). The temperature of the Simulated
%   Annealing is set to walk a large span of the parameter space to ensure
%   the global minimum is obtained. The minimisation function is based on
%   minimising the standard deviation of the sum of the squared residuals.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xpsStr:         MATLAB data-structure that contains the XPS data.
%   -   cTYPE:          1xN vector of the type of curve to use for the nth state.
%   -   iPESCurves:   	Nx8 array: the n'th peak parameters [BE,INT,FWHM,MR,LSE,LSI,LSW,ASY]
%   -   dPESCurves: 	Nx8 array: the n'th peak parameters uncertainties [BE,INT,FWHM,MR,LSE,LSI,LSW,ASY]
%   -   bTYPE:          string of the type of background to use for fitting. Default: "Poly" ("none", "Shir", "LinShir")
%   -   iPESBgrnd:      1x5 vector of the background parameters: [LHS,RHS,ORD,LAM,DEL,BGR]
%   -   dPESBgrnd:     	1x5 vector of the background parameter uncertainties: [LHS,RHS,ORD,LAM,DEL,BGR]
%   -   solve_type:     String of either "lsqnonlin", "simulannealbnd"
%
%   OUT:
%   -   fitStr:         MATLAB data-structure that contains all the fit parameters / variables / information

%% Default parameters
% - Default input parameters to use
if nargin < 8; solve_type = "simulannealbnd"; end
if nargin < 7; dPESBgrnd = [0.2, 0.2, 0, 0.5, 0.5, 0.05]; end
if nargin < 6; iPESBgrnd = [min(xpsStr.xdat(:))+0.2, max(xpsStr.xdat(:))-0.2, 1, 0.5, 0, 0]; end
if nargin < 5; bTYPE = "Poly"; end
if nargin < 4; dPESCurves = abs(0.10.*iPESCurves); dPESCurves(:,1) = 0.05; end
if isempty(solve_type); solve_type = "simulannealbnd"; end
if isempty(dPESBgrnd); dPESBgrnd = [0.2, 0.2, 0, 0.5, 0.5, 0.05]; end
if isempty(iPESBgrnd); iPESBgrnd = [min(xpsStr.xdat(:))+0.2, max(xpsStr.xdat(:))-0.2, 1, 0.5, 0, 0]; end
if isempty(bTYPE); bTYPE = "Poly"; end
if isempty(dPESCurves); dPESCurves = abs(0.10.*iPESCurves); dPESCurves(:,1) = 0.05; end
% -- Defining global variable that stores all relevant model and data variables
global XPSObj; 

%% - 1 - Loading in the XPS data to be fitted
XPSObj	= xpsStr;

%% - 2 - Appending the input arguments to the global variable
XPSObj.fit_args.solve_type  = solve_type;
XPSObj.fit_args.nSTATES   	= size(iPESCurves, 1);
XPSObj.fit_args.cTYPE    	= cTYPE;
XPSObj.fit_args.iPESCurves	= iPESCurves;
XPSObj.fit_args.dPESCurves	= dPESCurves;
XPSObj.fit_args.bTYPE    	= bTYPE;
XPSObj.fit_args.iPESBgrnd	= iPESBgrnd;
XPSObj.fit_args.dPESBgrnd	= dPESBgrnd;

%% - 3 - SIMULATED ANNEALING TO MINIMISE RESIDUALS AND DETERMINE BEST FIT 
% -- Defining the initial conditions and limits of the optimisation method
x0  = [iPESCurves(:); iPESBgrnd(:)];
lb = [iPESCurves(:)- dPESCurves(:); iPESBgrnd(:) - dPESBgrnd(:)];
ub = [iPESCurves(:)+ dPESCurves(:); iPESBgrnd(:) + dPESBgrnd(:)];
%% - 3.1 - LEAST SQUARES FITTING METHOD FOR NON-LINEAR FUNCTIONS
if solve_type == "lsqnonlin"
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions(@lsqnonlin,...
        'Algorithm', 'levenberg-marquardt', 'MaxFunctionEvaluations', 1e5,...
        'FunctionTolerance', 1e-8, 'OptimalityTolerance', 1e-8);
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    [params, ~, ~, ~, ~, ~] = lsqnonlin(@residual_function, x0, lb, ub, options);
%% - 3.2 - GLOBAL OPTIMISATION - SIMULATED ANNEALING TO MINIMISE RESIDUALS
elseif solve_type == "simulannealbnd"
    % -- Defining the optimisation options for the hybrid minimisation that takes place after convergence
    hybridopts = optimoptions( 'fmincon' , 'OptimalityTolerance' , 1e-8, 'MaxFunctionEvaluations', 1e5);
    % -- Defining the optimisation options for the simulated annealing method
    options = optimoptions(@simulannealbnd,...
        'InitialTemperature', 1e2, 'ReannealInterval', 1e2, 'TemperatureFcn', 'temperaturefast',...
        'MaxFunctionEvaluations', 1e5, 'FunctionTolerance', 1e-8, 'HybridFcn' , {'fmincon' , hybridopts});
    % to view annealing, add:
    % 'PlotFcn',{@saplotbestf,@saplottemperature,@saplotf,@saplotstopping}, 'MaxIterations', 50,
    % -- Run the simulated annealing method
    % --- (1) Initially optimise all parameters
    [params, ~, ~, ~] = simulannealbnd(@residual_function, x0, lb, ub, options);
    % --- (2) Keep all parameters constant, but optimise the peak heights, FWHM, MR only
    x0      = params; 
    lb2     = params; 
    ub2     = params;
    indx    = 1+XPSObj.fit_args.nSTATES:1:4*XPSObj.fit_args.nSTATES;
    lb2(indx) = lb(indx);
    ub2(indx) = ub(indx);
    % ----- Also allow constant background to vary
    lb2(end-5) = lb(end-5);
    ub2(end-5) = ub(end-5);
    [params, ~, ~, ~] = simulannealbnd(@residual_function, x0, lb2, ub2, options);
    % --- (3) Optimise all parameters again
    x0      = params;
    lb3     = lb;
    ub3     = ub;
    [params, ~, ~, ~] = simulannealbnd(@residual_function, x0, lb3, ub3, options);
end
%% - 4 - STORING THE FITS AS A MATLAB STRUCTURE
fitStr              = struct();
% -- Storing the initial arguments
fitStr.nSTATES      = size(iPESCurves, 1);
fitStr.cTYPE        = cTYPE;
fitStr.iPESCurves	= iPESCurves;
fitStr.dPESCurves	= dPESCurves;
fitStr.bTYPE        = bTYPE;
fitStr.iPESBgrnd	= iPESBgrnd;
fitStr.dPESBgrnd	= dPESBgrnd;
% -- Storing the final fit variables (PES curves)
fitStr.XX           = linspace(min(XPSObj.xdat(:)), max(XPSObj.xdat(:)), 1e3)';
fitStr.YY           = zeros(size(fitStr.XX));
for i = 1:fitStr.nSTATES
    % --- Best fit components on new domain
    fitStr.cPARAMS(i,:) = params(i:fitStr.nSTATES:end-6);
   	fitStr.cYY(:,i) = PESCurve(fitStr.XX, fitStr.cTYPE(i),...
        fitStr.cPARAMS(i,1), fitStr.cPARAMS(i,2), fitStr.cPARAMS(i,3),...
        fitStr.cPARAMS(i,4), fitStr.cPARAMS(i,5), fitStr.cPARAMS(i,6),...
        fitStr.cPARAMS(i,7), fitStr.cPARAMS(i,8));
    fitStr.YY = fitStr.YY + fitStr.cYY(:,i);
    % --- Storing each curve parameter
    fitStr.BE(i)    = fitStr.cPARAMS(i,1);
    fitStr.INT(i)   = fitStr.cPARAMS(i,2);
    fitStr.FWHM(i)  = fitStr.cPARAMS(i,3);
    fitStr.MR(i)    = fitStr.cPARAMS(i,4);
    fitStr.LSE(i)   = fitStr.cPARAMS(i,5);
    fitStr.LSI(i)   = fitStr.cPARAMS(i,6);
    fitStr.LSW(i)   = fitStr.cPARAMS(i,7);
    fitStr.ASY(i)   = fitStr.cPARAMS(i,8);
    % --- Storing the curve area information
    fitStr.AREA(i)	= trapz(fitStr.XX, fitStr.cYY(:,i));
end
% --- Storing the normalised curve area for quantification
fitStr.AREA0        = 100 * (fitStr.AREA ./ trapz(fitStr.XX, fitStr.YY));
% -- Storing the final fit variables (PES background)
fitStr.bPARAMS  	= params(end-5:end);
fitStr.ybgrnd       = PESBackground(XPSObj.xdat, XPSObj.int, fitStr.bTYPE,...
    fitStr.bPARAMS(1), fitStr.bPARAMS(2), fitStr.bPARAMS(3), fitStr.bPARAMS(4), fitStr.bPARAMS(5), fitStr.bPARAMS(6));
% -- Storing the initial data to be fitted
fitStr.xdat         = XPSObj.xdat;
fitStr.ydat_raw     = XPSObj.int;
fitStr.ydat         = XPSObj.int - fitStr.ybgrnd;
% -- Extract the residuals of the best curve fit
D   = fitStr.ydat_raw;
B   = fitStr.ybgrnd;
M 	= zeros(size(D));
for i = 1:fitStr.nSTATES
    cY(:,i) = PESCurve(XPSObj.xdat, fitStr.cTYPE(i),...
        fitStr.cPARAMS(i,1), fitStr.cPARAMS(i,2), fitStr.cPARAMS(i,3),...
        fitStr.cPARAMS(i,4), fitStr.cPARAMS(i,5), fitStr.cPARAMS(i,6),...
        fitStr.cPARAMS(i,7), fitStr.cPARAMS(i,8));
    M = M + cY(:,i);
end
fitStr.yresids      = (M - (D - B));
% -- Quantifying the quality of the fit with CHI-SQUARED (standard deviation of the residuals)
fitStr.CHISQ      	= std(fitStr.yresids);

end

%% DEFINING THE RESIDUAL MINIMISATION FUNCTION TO BE MINIMISED
function R = residual_function(x)
    global XPSObj;
    % - 1 - Extracting the DATA
    D = XPSObj.int;
    if size(D, 2) > 1; D = D'; end
    % - 2 - Extracting the MODEL
    M = model(x);
    if size(M, 2) > 1; M = M'; end
    % - 3 - Extracting the BACKGROUND
    B = background(x);
    if size(B, 2) > 1; B = B'; end
    % - 4 - Extracting the RESIDUALS
    R = std((M - (D - B)));
%     % - 5 - Plotting the lineshape & background as optimisation runs
%     figure(3633); cla reset; clf reset; hold on;
%     DB = D - B;
%     plot(XPSObj.xdat, DB, 'b-', 'linewidth', 2);
%     plot(XPSObj.xdat, M, 'k-', 'linewidth', 2);
%     gca_props();
%     xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
%     ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
%     axis([min(XPSObj.xdat(:)), max(XPSObj.xdat(:)), min(DB(:)), 1.2*max(DB(:))]);
end

%% DEFINING THE FUNCTION THAT DETERMINES THE TOTAL PES CURVE FIT
function M = model(x)
    % Calling the global variable
    global XPSObj;
    % - 1 - Extract the curve component for each state
    comp_int = {};
    for i = 1:XPSObj.fit_args.nSTATES
        % -- Extracting the arguments for the component curve
        pes_args    = x(i:XPSObj.fit_args.nSTATES:end-6);
        % -- Extracting the component intensities
        comp_int{i} = PESCurve(XPSObj.xdat, XPSObj.fit_args.cTYPE(i),...
            pes_args(1), pes_args(2), pes_args(3),...
            pes_args(4), pes_args(5), pes_args(6),...
            pes_args(7), pes_args(8));
    end
    % - 2 - Sum up all of the curve components to get the final model PES curve
    pes_int = zeros(size(XPSObj.xdat));
    for i = 1:XPSObj.fit_args.nSTATES
        pes_int = pes_int + comp_int{i};
    end
    M = pes_int;
end

%% DEFINING THE FUNCTION THAT DETERMINES THE TOTAL PES BACKGROUND
function B = background(x)
    % Calling the global variable
    global XPSObj;
    % - 1 - Determine the background to be used
    B = PESBackground(XPSObj.xdat, XPSObj.int, XPSObj.fit_args.bTYPE, x(end-5), x(end-4), x(end-3), x(end-2), x(end-1), x(end));
end