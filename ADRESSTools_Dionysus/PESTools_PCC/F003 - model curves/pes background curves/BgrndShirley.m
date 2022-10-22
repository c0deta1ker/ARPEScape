function [roi_xdat, roi_ydat, roi_bgrnd] = BgrndShirley(xdat, ydat, lhsVal, rhsVal, lhsWin, rhsWin, plot_result)
% [roi_xdat, roi_ydat, roi_bgrnd] = BgrndShirley(xdat, ydat, lhsVal, rhsVal, lhsWin, rhsWin, plot_result)
%   Function that determines the best fit to the background of PES (or 1D) 
%   data using a Shirley function. The output is a column vector containing
%   the background intensity, and the new domain and intensity of the data 
%   within the defined region of interest.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:           N×1 column vector of the input domain (binding energy for PES)
%   -   ydat:           N×1 column vector of the intensity range (intensity for PES)
%   -   lhsVal:         scalar of the LHS x-axis start position of the ROI.
%   -   rhsVal:         scalar of the RHS x-axis end position of the ROI.
%   -   lhsWin:         scalar of the LHS window size around the start position.
%   -   rhsWin:         scalar of the RHS window size around the end position.
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   roi_xdat:       M×1 column vector of the ROI x-domain (binding energy for PES).
%   -   roi_ydat:       M×1 column vector of the ROI y-values (intensity for PES).
%   -   roi_bgrnd:      M×1 column vector of the best fit Shirley background to the ROI data.

%% Default parameters
% - Setting the default conditions
if nargin < 3;          lhsVal = mean(xdat(:)) - abs(0.25*range(xdat(:))); end
if nargin < 4;          rhsVal = mean(xdat(:)) + abs(0.25*range(xdat(:))); end
if nargin < 5;          lhsWin = abs(0.05*range(xdat(:))); end
if nargin < 6;          rhsWin = abs(0.05*range(xdat(:))); end
if nargin < 7;          plot_result = 0; end
if isempty(lhsVal);     lhsVal = mean(xdat(:)) - abs(0.25*range(xdat(:))); end
if isempty(rhsVal);     rhsVal = mean(xdat(:)) + abs(0.25*range(xdat(:))); end
if isempty(lhsWin);     lhsWin = abs(0.05*range(xdat(:))); end
if isempty(rhsWin);     rhsWin = abs(0.05*range(xdat(:))); end
if isempty(plot_result); plot_result = 0; end
%% Validity checks on the input parameters
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

%% - 3 - Extracting the mean value of the step heights of the Shirley background
% -- Extracting a small range for averaging the intensity values
I1         	= mean(ydat(lhsWin_indx(1):lhsWin_indx(2)));
I2        	= mean(ydat(rhsWin_indx(1):rhsWin_indx(2)));

%% - 4 - Determination of the Shirley Background
% - Initialising the variables
Intensity = {};
ShirBgrnd = {};
% - Iterate 12 times to determine Shirley background
for N = 1:12
    % -- Create an initial guess of the background (use linear background)
    if N == 1
        [~, ~, ShirBgrnd{1}] = BgrndPoly(xdat, ydat, 1, lhsVal, rhsVal, lhsWin, rhsWin);
        Intensity{1} = roi_ydat - ShirBgrnd{1};
    % -- After initial guess, iterate to solution 
    else
        % --- File through all data points within the EB range
        for i = 1:nRows
            if i == 1
                A1	= 0;
                A2	= trapz(roi_xdat(i:end), Intensity{N-1}(i:end));
            elseif i == nRows
                A1	= trapz(roi_xdat(1:i), Intensity{N-1}(1:i));
                A2	= 0;
            else
                A1	= trapz(roi_xdat(1:i), Intensity{N-1}(1:i));
                A2	= trapz(roi_xdat(i:end), Intensity{N-1}(i:end));
            end
            ShirBgrnd{N}(i,1) = I2 + (I1 - I2) .* A2 ./ (A1 + A2);
        end
        Intensity{N} = roi_ydat - ShirBgrnd{N};
    end
end
% - Store the background intensity
roi_bgrnd = ShirBgrnd{end};
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
    gca_props(); title('BgrndShirley()'); 
    xlabel('$$ \bf  X $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Y $$', 'Interpreter', 'latex');
    legend({'Initial Data', 'ROI: Data', 'ROI: Background', 'ROI: Final'}, 'location', 'best', 'fontsize', 9);
    % -- Determining the best limits for the plot
    axLim_y = [roi_ydat; roi_ydat-roi_bgrnd];
    axis([min(xdat(:)), max(xdat(:)), min(axLim_y(:)), max(axLim_y(:))]);
end

end