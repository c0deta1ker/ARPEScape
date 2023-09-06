function crosshair_n(values, varargin)
% crosshair(values, varargin)
%   This function plots single cross-hair at the point defined by the user.
%   The properties of the line specifications can be defined by the user.
%
%   IN:
%   -   values:           1xN cell-array of the x- and y-positions of the cross-hairs.
%   -   varargin:         LineSpec arguments: LineSpec properties control the appearance and behavior of a LineSpec object, specified as a string scalar or character vector.
%
%   OUT:    (none)

%% Default parameters
if nargin < 2; varargin = {}; end
if isempty(varargin); varargin = {}; end
%% 1 - Plotting the vertical lines at the x-values
gca; hold on;
length(values)
for ii = 1:length(values)
    crosshair(values{ii}(1), values{ii}(2), varargin{:});
end
hold off