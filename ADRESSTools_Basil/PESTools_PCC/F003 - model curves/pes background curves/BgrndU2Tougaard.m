function [roi_xdat, roi_int, roi_bgrnd] = BgrndU2Tougaard(xdat, int, lhsLims, rhsLims, B, C)
% UNDER DEVELOPMENT!!!!!!!!!!!!!
% [roi_xdat, roi_int, roi_bgrnd] = BgrndU2Tougaard(xdat, int, lhsLims, rhsLims, B, C)
%   Function that determines the best fit to the background of PES data
%   using a U2 Tougaard function. The output in a colum vector containing
%   the background intensity, that is mapped onto the same domain as the
%   PES data. Calculated B and User defined C < 0 in U(x: B, C, 0, 0).
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       Nx1 column vector of the input domain (binding energy for XPS)
%   -   int:        Nx1 column vector of the intensity range (intensity for XPS)
%   -   lhsLims:   	value of the LHS binding energy position for the Tougaard backgorund
%   -   rhsLims: 	value of the RHS binding energy position for the Tougaard backgorund
%   -   C:          value of the constant C < 0 for U2 Tougaard background.
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

%% - 4 - Determination of the U2 Tougaard Background
% - Extracting the energy loss cross section F(x)
u2F             = @(x, B, C) (B.*x) ./ ((C + x.^2).^2);
F               = u2F(roi_xdat, B, C);
% - File through all data points within the EB range
T = [];
for i = 2:nRows
    t_xdat  = roi_xdat(1:i);
    t_int   = roi_int(1:i);
    t_F     = F(1:i);
    A       = trapz(t_xdat, t_F.*t_int);
    T(i) = A;
end

% - Store the background intensity
roi_bgrnd = T;
% -- Verifying outputs are all column vectors
if size(roi_xdat, 2) > 1;   roi_xdat = roi_xdat'; end
if size(roi_int, 2) > 1;    roi_int = roi_int'; end
if size(roi_bgrnd, 2) > 1;  roi_bgrnd = roi_bgrnd'; end


figure(); hold on;
plot(roi_xdat, roi_bgrnd);

end