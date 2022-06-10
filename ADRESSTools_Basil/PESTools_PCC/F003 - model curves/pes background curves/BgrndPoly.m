function [roi_xdat, roi_ydat, roi_bgrnd] = BgrndPoly(xdat, ydat, polyOrder, lhsVal, rhsVal, lhsWin, rhsWin, plot_result)
% [roi_xdat, roi_ydat, roi_bgrnd] = BgrndPoly(xdat, ydat, polyOrder, lhsVal, rhsVal, lhsWin, rhsWin, plot_result)
%   Function that determines the best fit to the background of PES (or 1D) 
%   data using a polynomial function. The output is a column vector containing
%   the background intensity, and the new domain and intensity of the data 
%   within the defined region of interest.
%   
%   REQ. FUNCTIONS: none
%   
%   IN:
%   -   xdat:           N×1 column vector of the input domain (binding energy for PES)
%   -   ydat:           N×1 column vector of the intensity range (intensity for PES)
%   -   polyOrder:      scalar, positive integer of the Polynomial Order.
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
if nargin < 3;          polyOrder = 1; end
if nargin < 4;          lhsVal = mean(xdat(:)) - abs(0.25*range(xdat(:))); end
if nargin < 5;          rhsVal = mean(xdat(:)) + abs(0.25*range(xdat(:))); end
if nargin < 6;          lhsWin = abs(0.05*range(xdat(:))); end
if nargin < 7;          rhsWin = abs(0.05*range(xdat(:))); end
if nargin < 8;          plot_result = 0; end
if isempty(polyOrder);  polyOrder = 1; end
if isempty(lhsVal);     lhsVal = mean(xdat(:)) - abs(0.25*range(xdat(:))); end
if isempty(rhsVal);     rhsVal = mean(xdat(:)) + abs(0.25*range(xdat(:))); end
if isempty(lhsWin);     lhsWin = abs(0.05*range(xdat(:))); end
if isempty(rhsWin);     rhsWin = abs(0.05*range(xdat(:))); end
if isempty(plot_result); plot_result = 0; end
%% Validity checks on the input parameters
% -- Making sure that the polynomial order is an integer value
polyOrder = floor(polyOrder(1));
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

%% - 3 - Determination of the Polynomial Background
% -- Extracting the data to be fitted with linear expressions
xdat_bgrnd      = [xdat(lhsWin_indx(1):lhsWin_indx(2));   xdat(rhsWin_indx(1):rhsWin_indx(2))];
ydat_bgrnd      = [ydat(lhsWin_indx(1):lhsWin_indx(2));   ydat(rhsWin_indx(1):rhsWin_indx(2))];
% -- Fitting the data to a polynomial of the correct order
fitType         = char("poly" + string(floor(polyOrder)));
if polyOrder == 0
    % -- For zero-order polynomial, only use a constant line
    roi_bgrnd       = mean(ydat_bgrnd(:)).*ones(size(roi_xdat));
else
    % -- Fitting the data to a polynomial of the correct order
    fitResult       = fit(xdat_bgrnd, ydat_bgrnd, fitType);
    % -- Extracting the background over a consistent domain
    roi_bgrnd       = fitResult(roi_xdat);
end
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
    gca_props(); title('BgrndPoly()'); 
    xlabel('$$ \bf  X $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Y $$', 'Interpreter', 'latex');
    legend({'Initial Data', 'ROI: Data', 'ROI: Background', 'ROI: Final'}, 'location', 'best', 'fontsize', 9);
    % -- Determining the best limits for the plot
    axLim_y = [roi_ydat; roi_ydat-roi_bgrnd];
    axis([min(xdat(:)), max(xdat(:)), min(axLim_y(:)), max(axLim_y(:))]);
end

end