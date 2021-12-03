function [roi_xdat, roi_int, roi_bgrnd] = BgrndShirley(xdat, int, lhsLims, rhsLims)
% [roi_xdat, roi_int, roi_bgrnd] = BgrndShirley(xdat, int, lhsLims, rhsLims)
%   Function that determines the best fit to the background of PES data
%   using a Shirley function. The output in a colum vector containing
%   the background intensity, that is mapped onto the same domain as the
%   PES data.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       Nx1 column vector of the input domain (binding energy for XPS)
%   -   int:        Nx1 column vector of the intensity range (intensity for XPS)
%   -   lhsLims:   	value of the LHS binding energy position for the Shirley backgorund
%   -   rhsLims: 	value of the RHS binding energy position for the Shirley backgorund
%
%   OUT:
%   -   roi_xdat:       Mx1 column vector of the ROI domain (binding energy for XPS).
%   -   roi_int:        Mx1 column vector of the ROI intensity (intensity for XPS).
%   -   roi_bgrnd:      Mx1 column vector of the best fit polynomial background to the ROI data.

%% Default parameters
% - Setting the default conditions
if nargin < 4; rhsLims = mean(xdat(:)) + 3; end
if nargin < 3; lhsLims = mean(xdat(:)) - 3; end
if isempty(rhsLims); rhsLims = mean(xdat(:)) + 3; end
if isempty(lhsLims); lhsLims = mean(xdat(:)) - 3; end
% - Sorting the limits of the background subtraction
lhsLims = sort(lhsLims);
rhsLims = sort(rhsLims);
% - Ensuring xdat and int are column vectors
if size(xdat, 2) > 1; xdat = xdat'; end
if size(int, 2) > 1; int = int'; end

%% - 1 - Extracting the START and END point of the background
% - Validity check on the LHS and RHS data points
% -- LHS validity
[~, lhsIndx]	= min(abs(xdat - lhsLims(1)));
% -- RHS validity
[~, rhsIndx]	= min(abs(xdat - rhsLims(1)));
% - ROI indexes
roi_indx    = [lhsIndx, rhsIndx];

%% - 2 - Cropping the data over the ROI whose background is extracted
roi_xdat    = xdat(roi_indx(1):roi_indx(2));
roi_int     = int(roi_indx(1):roi_indx(2));
nRows       = length(roi_xdat);

%% - 3 - Extracting the mean value of the step heights of the Shirley background
% -- Extracting a small range for averaging the intensity values
lhsAvgIndx      = lhsIndx + [-5, 5]; if lhsAvgIndx(1) < 1; lhsAvgIndx(1) = 1; end
rhsAvgIndx      = rhsIndx + [-5, 5]; if rhsAvgIndx(2) > length(xdat); rhsAvgIndx(2) = length(xdat); end
% -- Intensity at the START and END points
I1              = mean(int(lhsAvgIndx(1):lhsAvgIndx(2)));
I2              = mean(int(rhsAvgIndx(1):rhsAvgIndx(2)));

%% - 4 - Determination of the Shirley Background
% - Initialising the variables
Intensity = {};
ShirBgrnd = {};
% - Iterate 30 times to determine Shirley background
for N = 1:10
    % -- Create an initial guess of the background (use linear background)
    if N == 1
        [~, ~, ShirBgrnd{1}] = BgrndPoly(xdat, int, lhsLims, rhsLims, 1);
        Intensity{1} = roi_int - ShirBgrnd{1};
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
        Intensity{N} = roi_int - ShirBgrnd{N};
    end
end
% - Store the background intensity
roi_bgrnd = ShirBgrnd{end};
% -- Verifying outputs are all column vectors
if size(roi_xdat, 2) > 1;   roi_xdat = roi_xdat'; end
if size(roi_int, 2) > 1;    roi_int = roi_int'; end
if size(roi_bgrnd, 2) > 1;  roi_bgrnd = roi_bgrnd'; end

end