function [xVal, yVal] = find_peak_loc(xDat, yDat, xWin, type, plot_results)
% [xVal, yVal] = find_peak_loc(xDat, yDat, xWin, type, plot_results)
%   This function is used to determine both the maximum value and index of
%   a data set using a range of different methods; spline interpolation or
%   curve fitting a gaussian / voigt function. 
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   xdat:   	1D vector of the dataset domain.
%   -   ydat:    	1D vector of the dataset range.
%   -   xWin:     	Nx2 array that contains the x-windows where N peaks are expected to be found.
%   -   type:  		string of the type of method to use ("none","spline","gauss","smooth","xps_solver").
%   -   plot_results:  	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -   xVal:      	scalar or 1xN vector of the maxima x-values.
%   -   yVal:    	scalar or 1xN vector of the maxima y-values.

%% Default parameters
if nargin < 5; plot_results = 1; end 
if nargin < 4; type = "spline"; end 
if nargin < 3; xWin = []; type = "spline"; end 
if isempty(plot_results); plot_results = 1; end 
if isempty(xWin); xWin = []; end 
if isempty(type); type = "spline"; end
type = string(type);
pp = plot_props();

%% 1 - Initialising input parameters
if isempty(xWin)
    if min(xDat(:)) < 0;    minWin = 0.95*min(xDat(:)); 
    else;                   minWin = 1.05*min(xDat(:)); 
    end
    if max(xDat(:)) < 0;    maxWin = 1.05*max(xDat(:)); 
    else;                   maxWin = 0.95*max(xDat(:)); 
    end
    xWin	= [minWin, maxWin];
end
% % Linear interpolation of the data to ensure a dense grid
% ixDat   = linspace(min(xDat(:)), max(xDat(:)), 1e3);
% iyDat  	= interp1(xDat,yDat,ixDat, 'pchip');
% % Resetting data
% xDat = ixDat; yDat = iyDat;

%% 2 - Extracting the peak index and value
if type ~= "xps_solver"
    % -- Initialising arrays
    xVal = []; yVal = [];
    xx = {}; yy = {};
    % -- Extracting the indices over which to find the maximum
    for i = 1:size(xWin, 1)
        % -- Finding the indices
        [~, lbIndx]             = min(abs(xDat - xWin(i,1)));
        [~, ubIndx]             = min(abs(xDat - xWin(i,2)));
        Indx                    = [lbIndx, ubIndx];
        xDat_roi                = xDat(Indx(1):Indx(2));
        yDat_roi                = yDat(Indx(1):Indx(2));
        % -- Finding the maximum peak value from the data
        if type == "none"
            [yVal(i), maxInd]	= max(yDat_roi(:));
            xVal(i)           	= xDat_roi(maxInd);
        % -- Finding the maximum peak value from the smoothed data
        elseif type == "smooth"
            yDat_roi          	= Gaco1(yDat_roi, 2);
            [yVal(i), maxInd]	= max(yDat_roi(:));
            xVal(i)           	= xDat_roi(maxInd);
        % -- Finding the maximum peak value by fitting a spline
        elseif type == "spline"
            xx{i}             	= linspace(min(xDat_roi(:)), max(xDat_roi(:)), 1e3);
            yy{i}             	= spline(xDat_roi, yDat_roi, xx{i});
%             yy{i}             	= interp1(xDat_roi, yDat_roi, xx{i}, 'pchip');
            [yVal(i), maxInd]	= max(yy{i}(:));
            xVal(i)           	= xx{i}(maxInd);
        % -- Finding the maximum peak value by fitting a gaussian
        elseif type == "gauss"
            % Ensuring xdat is a column vector
            if size(xDat_roi, 2) >1; xDat_roi = xDat_roi'; end
            if size(yDat_roi, 2) >1; yDat_roi = yDat_roi'; end
            f       = fit(xDat_roi, yDat_roi, 'gauss1');
            xx{i}	= linspace(min(xDat_roi(:)), max(xDat_roi(:)), 1e3);
            yy{i}	= f(xx{i});
            [yVal(i), maxInd]	= max(yy{i}(:));
            xVal(i)           	= xx{i}(maxInd);
        end
    end
    % -- Initialising the figure
    if plot_results == 1
        fig = figure();
        fig.Position(3) = pp.fig4x4(1); 
        fig.Position(4) = pp.fig4x4(2);
        hold on;
        % Plotting the data
        plot(xDat, yDat, 'k.-', 'markersize', 12);
        % Plotting the outcome of the peak finder
        for i = 1:size(xWin, 1)
            plot(xVal(i), yVal(i), '.', 'markersize', 20, 'color', pp.col.fit{i},  'linewidth', 1.5);
            plot(xx{i}, yy{i}, '-', 'color', pp.col.fit{i});
            line([1 1].*xVal(i), [0, yVal(i)], 'Color', pp.col.fit{i}, 'LineWidth', 1.5, 'Linestyle', '-');
        end
        % - Formatting the figure
        gca_props(0); grid on;
        ylabel('$$ \bf  y $$', 'Interpreter', 'latex');
        xlabel('$$ \bf  x $$', 'Interpreter', 'latex');
        axis([min(xWin(:)), max(xWin(:)), 0.50.*min(yVal(:)), 1.25.*max(yVal(:))]);
        % Printing the peak finder values
        fprintf(("%i) x = %.2f eV, y = %.2f \n"), i, xVal(i), yVal(i));
    end
end

if type == "xps_solver"
    % (A) - Initialising the data structure
    xps_dat         = struct();
    xps_dat.xdat    = xDat;
    xps_dat.int     = yDat;
    xps_dat.raw_int = yDat;
    xps_dat.hv      = 0;
    n_curves        = size(xWin, 1);
    % - Normalising the peak height
    xps_dat.int = xps_dat.int - min(xps_dat.int(:));
    xps_dat.int = xps_dat.int ./ max(xps_dat.int(:));
    
    % (B) - DEFINING THE TYPES OF CURVES TO BE USED
    cTYPE   = strings(n_curves,1)+"sGLA";              % type of curve to use for fitting. Default: "sGLA" ("pGLA", "DS")
    % 1 :: DEFINING THE INITIAL CONDITIONS OF THE XPS COMPONENTS
    BE      = mean(xWin,2);           % scalar of the binding energy of PE curve. Each new row gives a new BE component. Make sure sizes are all consistent.
    INT     = ones(n_curves,1).*1;     % scalar of the peak intensity of PE curve.
    FWHM    = ones(n_curves,1).*0.05;  % scalar of the FWHM of the PE curve.
    MR      = ones(n_curves,1).*0.5;   % scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
    LSE     = zeros(n_curves,1);       % scalar of the binding energy of spin-orbit split PE curve.
    LSI     = zeros(n_curves,1);       % scalar of the branching ratio of spin-orbit split PE curve.
    LSW     = zeros(n_curves,1);       % scalar of the additional lorentzian width of spin-orbit split PE curve.
    ASY     = zeros(n_curves,1);       % scalar of the PE curve asymmetry parameter (usually for metallic systems).
    iparams = {}; iparams{1} = [BE, INT, FWHM, MR, LSE, LSI, LSW, ASY];
    % 2 :: DEFINING THE UNCERTAINTIES IN THE FIT PARAMETERS
    % -- Lower bounds
    iparams{2} = iparams{1}.*0; 
    iparams{2}(:,1) = xWin(:,1);
    iparams{2}(:,2) = 0; 
    iparams{2}(:,3) = 0; 
    iparams{2}(:,4) = 0;
    % -- Upper bounds
    iparams{3}      = iparams{1}.*0; 
    iparams{3}(:,1) = xWin(:,2);
    iparams{3}(:,2) = 10.0; 
    iparams{3}(:,3) = 2.0; 
    iparams{3}(:,4) = 1.0;  
    iparams{3}(:,8) = 1.0;  
    
    % (C) - Define the type of background to be used
    ibgrnd = {};
    % DEFINING THE TYPE OF BACKGROUND TO BE USED
    bTYPE   = "none";       % type of background to use for fitting. Default: "Poly" ("none", "Shir", "LinShir")
    % 1 :: DEFINING THE BACKGROUND PARAMETERS
    LHS     = min(xWin(:));        % scalar of the START point on the LHS for background
    RHS     = max(xWin(:));        % scalar of the END point on the RHS for background
    ORD     = 1;            % positive integer of the Polynomial Order for "Poly" background.
    LAM     = 0;          % scalar for "LinShir" mixing ratio: lambda = 0 is Pure Shirley; lambda = 1 is Pure Linear.
    DEL     = 0;            % scalar for the "LinShir" curve offset in binding energy.
    BGR     = 0;       % scalar for a constant background to be included in the fit
    ibgrnd{1} = [LHS, RHS, ORD, LAM, DEL, BGR];
    % 2 :: DEFINING THE UNCERTAINTY IN THE BACKGROUND PARAMETERS
    ibgrnd{2} = ibgrnd{1};  ibgrnd{2}(:,6) = ibgrnd{1}(:,6) - 0.05;
    ibgrnd{3} = ibgrnd{1}; ibgrnd{3}(:,6) = ibgrnd{1}(:,6) + 0.05;
    
    % (D) - Fit the data
    n_runs      = 1;
    solve_type  = "lsqcurvefit";
    xps_dat.fit = xps_solver(xps_dat, cTYPE, iparams, bTYPE, ibgrnd, solve_type);
    if plot_results == 1; view_xps_fit(xps_dat.fit); end
    
    % (E) - Assigning the data
    xVal = xps_dat.fit.BE;
    yVal = [];
    for i = 1:length(xVal)
        [~, indx] = min(abs(xps_dat.xdat - xVal(i)));
        yVal(i) = xps_dat.raw_int(indx);
    end
end
end