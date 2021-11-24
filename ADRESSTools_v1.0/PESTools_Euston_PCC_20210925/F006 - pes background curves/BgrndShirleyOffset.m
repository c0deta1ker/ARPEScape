function [roi_xdat, roi_int, roi_bgrnd] = BgrndShirleyOffset(xdat, int, lhsLims, rhsLims, lambda, delta)
% [roi_xdat, roi_int, roi_bgrnd] = BgrndShirleyOffset(xdat, int, lhsLims, rhsLims, lambda, delta)
%   Function that determines the best fit to the background of PES data
%   using an Offset Shirley function. This is a blend between a Shirley and
%   linear background and is defined by the parameters lambda and delta.
%   The output in a column vector containing the background intensity, 
%   that is mapped onto the same domain as the PES data.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       Nx1 column vector of the input domain (binding energy for XPS)
%   -   int:        Nx1 column vector of the intensity range (intensity for XPS)
%   -   lhsLims:   	value of the LHS binding energy position for the Shirley backgorund
%   -   rhsLims: 	value of the RHS binding energy position for the Shirley backgorund
%   -   lambda:   	scalar for mixing ratio: lambda = 0 is Pure Shirley; lambda = 1 is Pure Linear.
%   -   delta:      scalar for the Shirley curve offset in binding energy.
%
%   OUT:
%   -   roi_xdat:       Mx1 column vector of the ROI domain (binding energy for XPS).
%   -   roi_int:        Mx1 column vector of the ROI intensity (intensity for XPS).
%   -   roi_bgrnd:      Mx1 column vector of the best fit polynomial background to the ROI data.

%% Default parameters
% - Setting the default conditions
if nargin < 6; delta = 0; end
if nargin < 5; lambda = 0.5; end
if nargin < 4; rhsLims = mean(xdat(:)) + 3; end
if nargin < 3; lhsLims = mean(xdat(:)) - 3; end
if isempty(rhsLims); rhsLims = mean(xdat(:)) + 3; end
if isempty(lhsLims); lhsLims = mean(xdat(:)) - 3; end
if isempty(delta);  delta = 0; end
if isempty(lambda); lambda = 0.5; end 
% - Sorting the ROI
rhsLims = sort(rhsLims);
lhsLims = sort(lhsLims);
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

%% - 3 - Determination of the Offset Shirley Background
[~, ~, roi_bgrnd_lin]  	= BgrndPoly(xdat, int, min(roi_xdat(:)), max(roi_xdat(:)), 1);
[~, ~, roi_bgrnd_shrly]	= BgrndShirley(xdat - delta, int, min(roi_xdat(:)), max(roi_xdat(:)));
% - Determine the offset Shirley background
roi_bgrnd = roi_bgrnd_shrly.*(1-lambda) + lambda.*roi_bgrnd_lin;
% -- Verifying outputs are all column vectors
if size(roi_xdat, 2) > 1;   roi_xdat = roi_xdat'; end
if size(roi_int, 2) > 1;    roi_int = roi_int'; end
if size(roi_bgrnd, 2) > 1;  roi_bgrnd = roi_bgrnd'; end

end