function [roi_xdat, roi_int, roi_bgrnd] = BgrndPoly(xdat, int, lhsLims, rhsLims, Order)
% [roi_xdat, roi_int, roi_bgrnd] = BgrndPoly(xdat, int, lhsLims, rhsLims, Order)
%   Function that determines the best fit to the background of PES data
%   using a polynomial function. The output in a column vector containing
%   the background intensity, and the new domain and intensity of the XPS
%   data within the defined region of interest.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       Nx1 column vector of the input domain (binding energy for XPS)
%   -   int:        Nx1 column vector of the intensity range (intensity for XPS)
%   -   lhsLims:   	value of the LHS binding energy position for the Shirley backgorund
%   -   rhsLims: 	value of the RHS binding energy position for the Shirley backgorund
%   -   Order:     	positive integer of the Polynomial Order of background.
%
%   OUT:
%   -   roi_xdat:       Mx1 column vector of the ROI domain (binding energy for XPS).
%   -   roi_int:        Mx1 column vector of the ROI intensity (intensity for XPS).
%   -   roi_bgrnd:      Mx1 column vector of the best fit polynomial background to the ROI data.

%% Default parameters
% - Setting the default conditions
if nargin < 5; Order = 1; end
if nargin < 4; rhsLims = mean(xdat(:)) + 3; end
if nargin < 3; lhsLims = mean(xdat(:)) - 3; end
if isempty(rhsLims); rhsLims = mean(xdat(:)) + 3; end
if isempty(lhsLims); lhsLims = mean(xdat(:)) - 3; end 
if isempty(Order); Order = 1; end
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

%% - 3 - Determination of the Polynomial Background
% -- Extracting a small range for averaging the intensity values
lhsAvgIndx      = lhsIndx + [-5, 5]; if lhsAvgIndx(1) < 1; lhsAvgIndx(1) = 1; end
rhsAvgIndx      = rhsIndx + [-5, 5]; if rhsAvgIndx(2) > length(xdat); rhsAvgIndx(2) = length(xdat); end
% -- Extracting the data to be fitted with linear expressions
xdat_bgrnd      = [xdat(lhsAvgIndx(1):lhsAvgIndx(2));   xdat(rhsAvgIndx(1):rhsAvgIndx(2))];
int_bgrnd       = [int(lhsAvgIndx(1):lhsAvgIndx(2));    int(rhsAvgIndx(1):rhsAvgIndx(2))];
% -- Fitting the data to a polynomial of the correct order
fitType         = char("poly" + string(floor(Order)));
fitResult       = fit(xdat_bgrnd, int_bgrnd, fitType);
% -- Extracting the background over a consistent domain
roi_bgrnd 	= fitResult(roi_xdat);
% -- Verifying outputs are all column vectors
if size(roi_xdat, 2) > 1;   roi_xdat = roi_xdat'; end
if size(roi_int, 2) > 1;    roi_int = roi_int'; end
if size(roi_bgrnd, 2) > 1;  roi_bgrnd = roi_bgrnd'; end

end