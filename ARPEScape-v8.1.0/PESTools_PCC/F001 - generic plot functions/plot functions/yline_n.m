function yline_n(values, varargin)
% yline_n(values, varargin) 
%   This function plots single (or multiple) horizontal lines that cross the 
%   y-axis at the defined 'values'. The user can also choose the line style
%   of the lines.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   values:             1xN vector that contains all the y-values.
%   -   varargin:           LineSpec arguments: LineSpec properties control the appearance and behavior of a LineSpec object, specified as a string scalar or character vector.
%
%   OUT:    (none)

%% Default parameters
if nargin < 1; values = []; end
if nargin < 2; varargin = {}; end
if isempty(values); values = []; end
if isempty(varargin); varargin = {}; end
%% 1 - Plotting the vertical lines at the x-values
gca; hold on;
for ii = 1:length(values)
    yline(values(ii), varargin{:});
end
hold off