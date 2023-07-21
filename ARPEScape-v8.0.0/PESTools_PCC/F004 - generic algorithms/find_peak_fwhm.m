function [fwhm, fwhmLocs] = find_peak_fwhm(xdat, ydat, xWin, type, plot_result)
% [fwhm, fwhmLocs] = find_peak_fwhm(xdat, ydat, xWin, type, plot_result)
%   This function is used to determine both the fwhm of a data set using 
%   a range of different methods; spline interpolation or curve fitting a 
%   gaussian / voigt function. 
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   xdat:           [N×1] column vector of the input domain.
%   -   ydat:           [N×1] column vector of the output range.
%   -   xWin:           [N×2] array that contains the x-windows where N peaks are expected to be found.
%   -   type:           string of the type of method to use. Default: "spline" ("none","spline","gaco1","G","L","sGL","pGL","sGLA","pGLA").
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   fwhm:           scalar or [N×1] column vector of the fwhm of each peak defined.
%   -   fwhmLocs:   	[2N×2] column vector of the fwhm positions.

%% Default parameters
if nargin < 3; xWin = []; end
if nargin < 4; type = "spline"; end
if nargin < 5; plot_result = 0; end
if isempty(xWin); xWin = []; end
if isempty(type); type = "spline"; end
if isempty(plot_result); plot_result = 0; end
%% Validity checks on the input parameters
pp      = plot_props();
type    = string(type);
% -- Making the default window-size if it is undefined
if isempty(xWin)
    minWin  = mean(xdat(:)) - abs(0.25*range(xdat(:)));
    maxWin  = mean(xdat(:)) + abs(0.25*range(xdat(:)));
    xWin    = [minWin, maxWin];
end
% -- Making sure the window-size is an N×2 array
if size(xWin, 2) ~= 2; xWin = xWin'; end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end
if size(ydat, 2) > 1; ydat = ydat'; end

%% 1 - Extracting the peak position and values
xVal = []; yVal = [];
fwhm = []; fwhmLocs = {};
xdat_roi = {}; ydat_roi = {};
xx = {}; yy = {}; 
if length(xdat(:)) < 5e4; Nx = 5e4; else; Nx = length(xdat(:));end
for i = 1:size(xWin, 1)
    % -- Finding the indices and ROI
    [~, lbIndx]             = min(abs(xdat - xWin(i,1)));
    [~, ubIndx]             = min(abs(xdat - xWin(i,2)));
    Indx                    = [lbIndx, ubIndx];
    xdat_roi{i}             = xdat(Indx(1):Indx(2));
    ydat_roi{i}             = ydat(Indx(1):Indx(2));
    % -- Finding the maximum from the raw data
    if strcmpi(type,"none")
        xx{i}             	= xdat_roi{i};
        yy{i}             	= ydat_roi{i};
        [yVal(i), maxInd]	= max(yy{i}(:));
        xVal(i)           	= xx{i}(maxInd);
    % -- Finding the maximum peak value from the gaussian smoothed data
    elseif strcmpi(type,"gaco1")
        xx{i}             	= linspace(min(xdat_roi{i}(:)), max(xdat_roi{i}(:)), Nx);
        yy{i}             	= Gaco1(ydat_roi{i}, 2);
        yy{i}             	= interp1(xdat_roi{i}, yy{i}, xx{i});
        [yVal(i), maxInd]	= max(yy{i}(:));
        xVal(i)           	= xx{i}(maxInd);
    % -- Finding the maximum peak value by fitting a spline
    elseif strcmpi(type,"spline")
        xx{i}             	= linspace(min(xdat_roi{i}(:)), max(xdat_roi{i}(:)), Nx);
        yy{i}             	= spline(xdat_roi{i}, ydat_roi{i}, xx{i});
        [yVal(i), maxInd]	= max(yy{i}(:));
        xVal(i)           	= xx{i}(maxInd);
    % -- Finding the maximum peak value by fitting a gaussian
    elseif strcmpi(type,"G") || strcmpi(type,"gauss") || strcmpi(type,"gaussian")
        % --- Defining the fit function
        fit_func    = fittype(@(x0, peak, fwhm, x) Gauss(x, x0, peak, fwhm));
        x0_start    = mean(xdat_roi{i}(:));
        peak_start  = max(ydat_roi{i}(:));
        fwhm_start  = 0.5*range(xdat_roi{i}(:));
        % --- Executing the fitting operation
        [fit1,~,~]	= fit(xdat_roi{i},ydat_roi{i},fit_func,'start',[x0_start, peak_start, fwhm_start]);
        params    	= coeffvalues(fit1);
        xx{i}   = linspace(min(xdat_roi{i}(:)), max(xdat_roi{i}(:)), 1e3);
        yy{i}	= Gauss(xx{i}, params(1), params(2), params(3));
        [yVal(i), maxInd]	= max(yy{i}(:));
        xVal(i)           	= xx{i}(maxInd);
    % -- Finding the maximum peak value by fitting a lorentzian
    elseif strcmpi(type,"L") || strcmpi(type,"lorz") || strcmpi(type,"lorentzian")
        % --- Defining the fit function
        fit_func    = fittype(@(x0, peak, fwhm, x) Lorentzian(x, x0, peak, fwhm));
        x0_start    = mean(xdat_roi{i}(:));
        peak_start  = max(ydat_roi{i}(:));
        fwhm_start  = 0.5*range(xdat_roi{i}(:));
        % --- Executing the fitting operation
        [fit1,~,~]	= fit(xdat_roi{i},ydat_roi{i},fit_func,'start',[x0_start, peak_start, fwhm_start]);
        params    	= coeffvalues(fit1);
        xx{i}   = linspace(min(xdat_roi{i}(:)), max(xdat_roi{i}(:)), Nx);
        yy{i}	= Lorentzian(xx{i}, params(1), params(2), params(3));
        [yVal(i), maxInd]	= max(yy{i}(:));
        xVal(i)           	= xx{i}(maxInd);
    % -- Finding the maximum peak value by fitting a summed voigt
    elseif strcmpi(type,"sGL")
        % --- Defining the fit function
        fit_func    = fittype(@(x0, peak, fwhm, mr, x) sGL(x, x0, peak, fwhm, mr));
        x0_start    = mean(xdat_roi{i}(:));
        peak_start  = max(ydat_roi{i}(:));
        fwhm_start  = 0.5*range(xdat_roi{i}(:));
        mr_start    = 0.5;
        % --- Executing the fitting operation
        [fit1,~,~]	= fit(xdat_roi{i},ydat_roi{i},fit_func,'start',[x0_start, peak_start, fwhm_start, mr_start]);
        params    	= coeffvalues(fit1);
        xx{i}   = linspace(min(xdat_roi{i}(:)), max(xdat_roi{i}(:)), Nx);
        yy{i}	= sGL(xx{i}, params(1), params(2), params(3), params(4));
        [yVal(i), maxInd]	= max(yy{i}(:));
        xVal(i)           	= xx{i}(maxInd);
    % -- Finding the maximum peak value by fitting a multiplied voigt
    elseif strcmpi(type,"pGL")
        % --- Defining the fit function
        fit_func    = fittype(@(x0, peak, fwhm, mr, x) pGL(x, x0, peak, fwhm, mr));
        x0_start    = mean(xdat_roi{i}(:));
        peak_start  = max(ydat_roi{i}(:));
        fwhm_start  = 0.5*range(xdat_roi{i}(:));
        mr_start    = 0.5;
        % --- Executing the fitting operation
        [fit1,~,~]	= fit(xdat_roi{i},ydat_roi{i},fit_func,'start',[x0_start, peak_start, fwhm_start, mr_start]);
        params    	= coeffvalues(fit1);
        xx{i}   = linspace(min(xdat_roi{i}(:)), max(xdat_roi{i}(:)), Nx);
        yy{i}	= pGL(xx{i}, params(1), params(2), params(3), params(4));
        [yVal(i), maxInd]	= max(yy{i}(:));
        xVal(i)           	= xx{i}(maxInd);
    % -- Finding the maximum peak value by fitting a summed, asymmetric voigt
    elseif strcmpi(type,"sGLA")
        % --- Defining the fit function
        fit_func    = fittype(@(x0, peak, fwhm, mr, asym, x) sGLA(x, x0, peak, fwhm, mr, asym));
        x0_start    = mean(xdat_roi{i}(:));
        peak_start  = max(ydat_roi{i}(:));
        fwhm_start  = 0.5*range(xdat_roi{i}(:));
        mr_start    = 0.5;
        asym_start  = 0.2;
        % --- Executing the fitting operation
        [fit1,~,~]	= fit(xdat_roi{i},ydat_roi{i},fit_func,'start',[x0_start, peak_start, fwhm_start, mr_start, asym_start]);
        params    	= coeffvalues(fit1);
        xx{i}   = linspace(min(xdat_roi{i}(:)), max(xdat_roi{i}(:)), Nx);
        yy{i}	= sGLA(xx{i}, params(1), params(2), params(3), params(4), params(5));
        [yVal(i), maxInd]	= max(yy{i}(:));
        xVal(i)           	= xx{i}(maxInd);
    % -- Finding the maximum peak value by fitting a multiplied, asymmetric voigt
    elseif strcmpi(type,"pGLA")
        % --- Defining the fit function
        fit_func    = fittype(@(x0, peak, fwhm, mr, asym, x) pGLA(x, x0, peak, fwhm, mr, asym));
        x0_start    = mean(xdat_roi{i}(:));
        peak_start  = max(ydat_roi{i}(:));
        fwhm_start  = 0.5*range(xdat_roi{i}(:));
        mr_start    = 0.5;
        asym_start  = 0.2;
        % --- Executing the fitting operation
        [fit1,~,~]	= fit(xdat_roi{i},ydat_roi{i},fit_func,'start',[x0_start, peak_start, fwhm_start, mr_start, asym_start]);
        params    	= coeffvalues(fit1);
        xx{i}   = linspace(min(xdat_roi{i}(:)), max(xdat_roi{i}(:)), Nx);
        yy{i}	= pGLA(xx{i}, params(1), params(2), params(3), params(4), params(5));
        [yVal(i), maxInd]	= max(yy{i}(:));
        xVal(i)           	= xx{i}(maxInd);
    end
    % -- Extracting the FWHM of the data
    yHalf(i)            = 0.5*yVal(i);                          % Find the half max value
    index1              = find(yy{i}(:) >= yHalf(i), 1, 'first');  % Find where the data first drops below half the max
    index2              = find(yy{i}(:) >= yHalf(i), 1, 'last');   % Find where the data last rises above half the max
    fwhm(i)             = abs(xx{i}(index2) - xx{i}(index1));   % Extracting FWHM
    fwhmLocs{i,1}    	= [xx{i}(index1), yy{i}(index1); xx{i}(index2), yy{i}(index2)];
end
fwhmLocs                = cell2mat(fwhmLocs);
%% 2 - Ensuring the peak positions and location are column vectors
if size(fwhm, 2) > 1; fwhm = fwhm'; end

%% -- For Debugging
if plot_result == 1
    fig = figure();
    fig.Position(3) = pp.fig4x4(1); 
    fig.Position(4) = pp.fig4x4(2);
    hold on;
    % Plotting the data
    plot(xdat, ydat, 'b.-', 'linewidth', 0.5);
    % Plotting the outcome of the peak finder
    for i = 1:size(xWin, 1)
        plot(xdat_roi{i}, ydat_roi{i}, 'r-','linewidth', 2.0);
        plot(xx{i}, yy{i}, 'k-', 'linewidth', 2.5, 'Color', pp.col.fit{i});
        % -- Plotting peak position
        plot(xVal(i), yVal(i), '.', 'markersize', 20, 'color', pp.col.fit{i},  'linewidth', 2.5);
        line([1 1].*xVal(i), [0, yVal(i)], 'Color', pp.col.fit{i}, 'LineWidth', 2.5, 'Linestyle', '-');
        % -- Plotting LHS point
        plot(fwhmLocs(2*i-1,1), fwhmLocs(2*i-1,2), '.', 'markersize', 20, 'color', pp.col.fit{i},  'linewidth', 2.5);
        line([1 1].*fwhmLocs(2*i-1,1), [0, yVal(i)], 'Color', pp.col.fit{i}, 'LineWidth', 1.0, 'Linestyle', '--');
        % -- Plotting RHS point
        plot(fwhmLocs(2*i,1), fwhmLocs(2*i,2), '.', 'markersize', 20, 'color', pp.col.fit{i},  'linewidth', 2.5);
        line([1 1].*fwhmLocs(2*i,1), [0, yVal(i)], 'Color', pp.col.fit{i}, 'LineWidth', 1.0, 'Linestyle', '--');
        % -- Printing useful information
        fprintf(("%i) x = %.2f eV, y = %.2f, fwhm = %.2f \n"), i, xVal(i), yVal(i), fwhm(i));
    end
    % - Formatting the figure
    gca_props(0); grid on; title('find_peak_fwhm()', 'interpreter', 'none'); 
    xlabel('$$ \bf  X $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Y $$', 'Interpreter', 'latex');
    axis([mean(xWin(:))-0.75*range(xWin(:)), mean(xWin(:))+0.75*range(xWin(:)), min(ydat(:)), 1.10.*max(yVal(:))]);
end

end