function [xVal, yVal] = find_peak_loc(xdat, ydat, xWin, type, plot_result)
% [xVal, yVal] = find_peak_loc(xdat, ydat, xWin, type, plot_results)
%   This function is used to determine both the maximum value and position of
%   a data set using a range of different methods; spline interpolation or
%   curve fitting a gaussian / voigt function. 
%
%   IN:
%   -   xdat:           [N×1] column vector of the input domain.
%   -   ydat:           [N×1] column vector of the output range.
%   -   xWin:           [N×2] array that contains the x-windows where N peaks are expected to be found.
%   -   type:           string of the type of method to use. Default: "spline" ("none","spline","gaco1","G","L","sGL","pGL","sGLA","pGLA").
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   xVal:           scalar or [N×1] column vector of the maxima x-values.
%   -   yVal:           scalar or [N×1] column vector of the maxima y-values.

%% Default parameters
if nargin < 3; xWin = []; end
if nargin < 4; type = ""; end
if nargin < 5; plot_result = 0; end
if isempty(xWin); xWin = []; end
if isempty(type); type = ""; end
if isempty(plot_result); plot_result = 0; end
%% Validity checks on the input parameters
pp      = plot_props();
type    = string(type);
% -- Making the default window-size if it is undefined
if isempty(xWin)
    minWin  = mean(xdat(:)) - abs(0.475*range(xdat(:)));
    maxWin  = mean(xdat(:)) + abs(0.475*range(xdat(:)));
    xWin    = [minWin, maxWin];
end
% -- Making sure the window-size is an N×2 array
if size(xWin, 2) ~= 2; xWin = xWin'; end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end
if size(ydat, 2) > 1; ydat = ydat'; end

%% 1 - Extracting the peak position and values
Npoints = 5e3;
xVal = []; yVal = [];
xdat_roi = {}; ydat_roi = {};
xx = {}; yy = {};
for i = 1:size(xWin, 1)
    % -- Finding the indices and ROI
    [~, lbIndx]             = min(abs(xdat - xWin(i,1)));
    [~, ubIndx]             = min(abs(xdat - xWin(i,2)));
    Indx                    = [lbIndx, ubIndx];
    xdat_roi{i}             = xdat(Indx(1):Indx(2));
    ydat_roi{i}             = ydat(Indx(1):Indx(2));
    if length(xdat_roi{i}(:)) < Npoints;  xx{i} = linspace(min(xdat_roi{i}(:)), max(xdat_roi{i}(:)), Npoints);
    else; xx{i} = linspace(min(xdat_roi{i}(:)), max(xdat_roi{i}(:)), length(xdat_roi{i}(:)));
    end
    % -- Finding the maximum peak value from the raw data
    if strcmpi(type,"none") || strcmpi(type,"")
        xx{i}             	= xdat_roi{i};
        yy{i}             	= ydat_roi{i};
        [yVal(i), maxInd]	= max(yy{i}(:));
        xVal(i)           	= xx{i}(maxInd);
    % -- Finding the maximum peak value from the gaussian smoothed data
    elseif strcmpi(type,"gaco1")
        ceil(0.05*length(ydat_roi{i}))
        yy{i}             	= Gaco1(ydat_roi{i}, ceil(0.02*length(ydat_roi{i})));
        yy{i}             	= interp1(xdat_roi{i}, yy{i}, xx{i});
        [yVal(i), maxInd]	= max(yy{i}(:));
        xVal(i)           	= xx{i}(maxInd);
    % -- Finding the maximum peak value by fitting a spline
    elseif strcmpi(type,"spline")
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
        yy{i}	= pGLA(xx{i}, params(1), params(2), params(3), params(4), params(5));
        [yVal(i), maxInd]	= max(yy{i}(:));
        xVal(i)           	= xx{i}(maxInd);
    end
end
%% 2 - Ensuring the peak positions and location are column vectors
if size(xVal, 2) > 1; xVal = xVal'; end
if size(yVal, 2) > 1; yVal = yVal'; end

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
        plot(xVal(i), yVal(i), '.', 'markersize', 20, 'color', pp.col.fit{i},  'linewidth', 2.5);
        line([1 1].*xVal(i), [0, yVal(i)], 'Color', pp.col.fit{i}, 'LineWidth', 2.5, 'Linestyle', '-');
        fprintf(("%i) x = %.2f, y = %.2f \n"), i, xVal(i), yVal(i));
    end
    % - Formatting the figure
    gca_props(0); title('find_peak_loc()', 'interpreter', 'none'); 
%     grid on;
    xlabel(' X ', 'fontweight', 'bold');
    ylabel(' Y ', 'fontweight', 'bold');
    if min(ydat(:)) ~= max(yVal(:))
        axis([mean(xWin(:))-0.75*(range(xWin(:))), mean(xWin(:))+0.75*(range(xWin(:))),...
            min(cell2mat(ydat_roi(:))), 1.10.*max(yVal(:))]);
    end
end

end