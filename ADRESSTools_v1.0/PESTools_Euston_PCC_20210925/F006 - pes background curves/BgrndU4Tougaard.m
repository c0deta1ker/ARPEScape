function [roi_xdat, roi_int, F] = BgrndU4Tougaard(xdat, int, lhsLims, rhsLims, B, C, t0, D)
% UNDER DEVELOPMENT
% [roi_xdat, roi_int, roi_bgrnd] = BgrndU4Tougaard(xdat, int, lhsLims, rhsLims, B, C, t0, D)
%   Function that determines the best fit to the background of PES data
%   using a U4 Tougaard function. The output in a colum vector containing
%   the background intensity, that is mapped onto the same domain as the
%   PES data. For the U4 Tougaard background, the parameters B, C, D and t0
%   are user defined in U(x: B, C, D, t0).
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       Nx1 column vector of the input domain (binding energy for XPS)
%   -   int:        Nx1 column vector of the intensity range (intensity for XPS)
%   -   lhsLims:   	value of the LHS binding energy position for the Tougaard backgorund
%   -   rhsLims: 	value of the RHS binding energy position for the Tougaard backgorund
%   -   B:          user defined B for U4 Tougaard background.
%   -   C:          user defined C for U4 Tougaard background.
%   -   t0:         user defined t0 for U4 Tougaard background.
%   -   D:          user defined D for U4 Tougaard background.
%
%   OUT:
%   -   roi_xdat:       Mx1 column vector of the ROI domain (binding energy for XPS).
%   -   roi_int:        Mx1 column vector of the ROI intensity (intensity for XPS).
%   -   roi_bgrnd:      Mx1 column vector of the best fit polynomial background to the ROI data.

%% Default parameters
% - Setting the default conditions
if nargin < 8; D = 0; end
if nargin < 7; t0 = 0; end
if nargin < 6; C = 0; end
if nargin < 5; B = 0; end
if nargin < 4; rhsLims = mean(xdat(:)) + 3; end
if nargin < 3; lhsLims = mean(xdat(:)) - 3; end
if isempty(D); D = 0; end
if isempty(t0); t0 = 0; end
if isempty(C); C = 0; end
if isempty(B); B = 0; end
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

%% - 4 - Determination of the U4 Tougaard Background
% - Extracting the energy loss cross section F(x)
u4F             = @(x, B, C, D) (B.*x) ./ ((C - x.^2).^2 + (D.*x.^2));
F               = u4F(roi_xdat, B, C, D);
[~, t0Indx]     = min(abs(roi_xdat - t0));
F(t0Indx:end)   = 0;        % F(x) = 0 when x >= t0
% - Determination of the Tougaard background
T               = trapz(roi_xdat, F.*roi_int);
% - Store the background intensity
roi_bgrnd = [0, T];
% -- Verifying outputs are all column vectors
if size(roi_xdat, 2) > 1;   roi_xdat = roi_xdat'; end
if size(roi_int, 2) > 1;    roi_int = roi_int'; end
if size(roi_bgrnd, 2) > 1;  roi_bgrnd = roi_bgrnd'; end

end
