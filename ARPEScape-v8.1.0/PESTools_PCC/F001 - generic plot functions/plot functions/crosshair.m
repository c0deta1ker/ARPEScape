function crosshair(xval, yval, varargin)
% crosshair(xval, yval, varargin)
%   This function plots single cross-hair at the point defined by the user.
%   The properties of the line specifications can be defined by the user.
%
%   IN:
%   -   xval:             single value of the x-position of the cross-hair.
%   -   yval:             single value of the y-position of the cross-hair.
%   -   varargin:         LineSpec arguments: LineSpec properties control the appearance and behavior of a LineSpec object, specified as a string scalar or character vector.
%
%   OUT:    (none)

%% Default parameters
if nargin < 3; varargin = {}; end
if isempty(varargin); varargin = {}; end
%% 1 - Plotting the vertical lines at the x-values
gca; hold on;
xline(xval, varargin{:});
yline(yval, varargin{:});
hold off