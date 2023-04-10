function [roi_xdat, roi_ydat, roi_bgrnd] = BgrndStepFDDGpL(xdat, ydat, fdd_ef, fdd_T, fdd_fwhm, lhsVal, rhsVal, lhsWin, rhsWin, plot_result)
% [roi_xdat, roi_ydat, roi_bgrnd] = BgrndStepFDDGpL(xdat, ydat, fdd_ef, fdd_T, fdd_fwhm, lhsVal, rhsVal, lhsWin, rhsWin, plot_result)
%   Function that determines the best fit to the background of PES (or 1D) 
%   data using a single FDD step. The output is a column vector containing
%   the background intensity, and the new domain and intensity of the data 
%   within the defined region of interest. Here, a linear background is
%   MULTIPLIED to the FDD function, so the background is fitted to the
%   the gradient, intercept and constant.
%   
%   REQ. FUNCTIONS: none
%   
%   IN:
%   -   xdat:           N×1 column vector of the input domain (binding energy for PES)
%   -   ydat:           N×1 column vector of the intensity range (intensity for PES)
%   -   fdd_ef:         scalar of the Fermi-level position (eV).
%   -   fdd_T:          scalar of the temperature (K).
%   -   fdd_fwhm:       scalar of the full-width at half-maximum (FWHM) of the Gaussian broadening term (eV).
%   -   lhsVal:         scalar of the LHS x-axis start position of the ROI.
%   -   rhsVal:         scalar of the RHS x-axis end position of the ROI.
%   -   lhsWin:         scalar of the LHS window size around the start position.
%   -   rhsWin:         scalar of the RHS window size around the end position.
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   roi_xdat:       M×1 column vector of the ROI x-domain (binding energy for PES).
%   -   roi_ydat:       M×1 column vector of the ROI y-values (intensity for PES).
%   -   roi_bgrnd:      M×1 column vector of the best fit polynomial background to the ROI data.

%% Default parameters
% - Setting the default conditions 
if nargin < 3;          fdd_ef      = 0.0; end
if nargin < 4;          fdd_T       = 12; end
if nargin < 5;          fdd_fwhm    = 0.1; end
if nargin < 6;          lhsVal = mean(xdat(:)) - abs(0.25*range(xdat(:))); end
if nargin < 7;          rhsVal = mean(xdat(:)) + abs(0.25*range(xdat(:))); end
if nargin < 8;          lhsWin = abs(0.02*range(xdat(:))); end
if nargin < 9;          rhsWin = abs(0.02*range(xdat(:))); end
if nargin < 10;         plot_result = 0; end
if isempty(fdd_ef);   	fdd_ef      = 0.0; end
if isempty(fdd_T);    	fdd_T       = 12; end
if isempty(fdd_fwhm); 	fdd_fwhm    = 0.1; end
if isempty(lhsVal);     lhsVal = mean(xdat(:)) - abs(0.25*range(xdat(:))); end
if isempty(rhsVal);     rhsVal = mean(xdat(:)) + abs(0.25*range(xdat(:))); end
if isempty(lhsWin);     lhsWin = abs(0.02*range(xdat(:))); end
if isempty(rhsWin);     rhsWin = abs(0.02*range(xdat(:))); end
if isempty(plot_result); plot_result = 0; end
%% Validity checks on the input parameters
% -- Validity check on the FDD parameters
if fdd_T < 0; fdd_T = 0; end
if fdd_fwhm < 0; fdd_fwhm = 0; end
% -- If a vector / array of the LHS and RHS are entered, only select the first value
if length(lhsVal) > 1
    disp('Warning: lhsVal is defined as a vector, only choosing the first element.');
    lhsVal = lhsVal(1);
end
if length(rhsVal) > 1
    disp('Warning: rhsVal is defined as a vector, only choosing the first element.');
    rhsVal = rhsVal(1);
end
% -- If the LHS value is smaller than the RHS value, switch the assignment
if lhsVal > rhsVal
    disp('Warning: lhsVal > rhsVal, values have been switched.');
    A = lhsVal; Awin = lhsWin;
    B = rhsVal; Bwin = rhsWin;
    lhsVal = B; lhsVal = Bwin;
    rhsVal = A; rhsVal = Awin;
% -- If the LHS and RHS values are identical, no window can be defined
elseif lhsVal == rhsVal
    error('Error: lhsVal == rhsVal, make sure you define a finite window size for the ROI.');
end
% -- If a vector / array of the LHS and RHS windows are entered, only select the first value
if length(lhsWin) > 1
    disp('Warning: lhsWin is defined as a vector, only choosing the first element.');
    lhsWin = lhsWin(1);
end
if length(rhsWin) > 1
    disp('Warning: rhsWin is defined as a vector, only choosing the first element.');
    rhsWin = rhsWin(1);
end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end
if size(ydat, 2) > 1; ydat = ydat'; end

%% - 1 - Extracting the index START and END point of the background
% -- LHS index points
[~, lhsIndx_mu]     = min(abs(xdat - lhsVal));
[~, lhsIndx_lb]     = min(abs(xdat - lhsVal - lhsWin));
[~, lhsIndx_ub]     = min(abs(xdat - lhsVal + lhsWin));
% -- RHS index points
[~, rhsIndx_mu]     = min(abs(xdat - rhsVal));
[~, rhsIndx_lb]     = min(abs(xdat - rhsVal - rhsWin));
[~, rhsIndx_ub]     = min(abs(xdat - rhsVal + rhsWin));
% - ROI indexes
roi_indx            = sort([lhsIndx_mu, rhsIndx_mu]);
lhsWin_indx      	= sort([lhsIndx_lb, lhsIndx_ub]);
rhsWin_indx     	= sort([rhsIndx_lb, rhsIndx_ub]);

%% - 2 - Cropping the data over the ROI whose background is extracted
roi_xdat    = xdat(roi_indx(1):roi_indx(2));
roi_ydat    = ydat(roi_indx(1):roi_indx(2));
nRows       = length(roi_xdat);

%% - 3 - Determination & fitting of the FDD Background
% -- Extracting the data to be fitted with linear expressions
xdat_bgrnd      = [xdat(lhsWin_indx(1):lhsWin_indx(2));   xdat(rhsWin_indx(1):rhsWin_indx(2))];
ydat_bgrnd      = [ydat(lhsWin_indx(1):lhsWin_indx(2));   ydat(rhsWin_indx(1):rhsWin_indx(2))];
% -- Fitting the data to a polynomial of the correct order
fermi_fit = fittype(@(a, b, c, x) FDDGpL(x, fdd_ef, fdd_T, fdd_fwhm, a, b, c));
a_start = 0.1;      % lin_grad:         scalar of the gradient of the linear background.
b_start = 1.0;      % lin_offset:       scalar of the y-intercept of the linear background.
c_start = 0.1;      % const_bgrnd:  	scalar of the constant background y-offset value.
% Executing the fitting operation
[fit1,~,~]          = fit(xdat_bgrnd,ydat_bgrnd,fermi_fit,'start',[a_start, b_start, c_start]);
fit1para            = coeffvalues(fit1);
roi_bgrnd           = FDDGpL(roi_xdat, fdd_ef, fdd_T, fdd_fwhm, fit1para(1), fit1para(2), fit1para(3));
% -- Verifying outputs are all column vectors
if size(roi_xdat, 2) > 1;   roi_xdat = roi_xdat'; end
if size(roi_ydat, 2) > 1;   roi_ydat = roi_ydat'; end
if size(roi_bgrnd, 2) > 1;  roi_bgrnd = roi_bgrnd'; end

%% -- For Debugging
if plot_result == 1
    % - Initialising the figure object
    figure(); hold on;
    % -- Plotting the ROI, LHS and RHS analysis windows
    ROI     = [xdat(roi_indx(1)), xdat(roi_indx(2))];
    LHS     = [xdat(lhsWin_indx(1)), xdat(lhsWin_indx(2))];
    RHS     = [xdat(rhsWin_indx(1)), xdat(rhsWin_indx(2))];
    hLHS    = patch([LHS(1), LHS(1), LHS(2), LHS(2), LHS(1)], [-1, 1, 1, -1, -1].*1e4, [0.6 0.6 0.2], 'facealpha', 0.4, 'edgecolor', 'none');
    hLHS.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hRHS    = patch([RHS(1), RHS(1), RHS(2), RHS(2), RHS(1)], [-1, 1, 1, -1, -1].*1e4, [0.6 0.6 0.2], 'facealpha', 0.4, 'edgecolor', 'none');
    hRHS.Annotation.LegendInformation.IconDisplayStyle = 'off';
    hMID    = patch([ROI(1), ROI(1), ROI(2), ROI(2), ROI(1)], [-1, 1, 1, -1, -1].*1e4, [0.8 0.9 0.8], 'facealpha', 0.5, 'edgecolor', 'none');
    hMID.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % -- Plotting a vertical line to show the ROI
    a = line([ROI(1) ROI(1)], [-1e5, 1e5], 'Color', [0 0 0], 'LineWidth', 1, 'Linestyle', '-');
    a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    b = line([ROI(2) ROI(2)], [-1e5, 1e5], 'Color', [0 0 0], 'LineWidth', 1, 'Linestyle', '-');
    b.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % -- Plotting the 1D data
    plot(xdat, ydat, 'b-', 'linewidth', 0.5);
    plot(roi_xdat, roi_ydat, 'b-', 'linewidth', 2);
    plot(roi_xdat, roi_bgrnd, 'r-', 'linewidth', 2);
    plot(roi_xdat, roi_ydat-roi_bgrnd, 'k-', 'linewidth', 2);
    gca_props(); title('BgrndStepFDDGpL()'); 
    xlabel('$$ \bf  X $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Y $$', 'Interpreter', 'latex');
    legend({'Initial Data', 'ROI: Data', 'ROI: Background', 'ROI: Final'}, 'location', 'best', 'fontsize', 9);
    % -- Determining the best limits for the plot
    axLim_y = [roi_ydat; roi_ydat-roi_bgrnd];
    axis([min(xdat(:)), max(xdat(:)), min(axLim_y(:)), max(axLim_y(:))]);
end

end