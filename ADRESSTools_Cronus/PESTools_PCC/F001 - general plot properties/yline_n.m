function yline_n(values,line_style) 
% yline_n(values,line_style) 
%   This function plots single (or multiple) horizontal lines that cross the 
%   y-axis at the defined 'values'. The user can also choose the line style
%   of the lines.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   values:     1 x N vector that contains all the y-values.
%   -   line_style:	char of the line style to be used for the plot.
%
%   OUT:    (none)

%% Default parameters
if nargin < 1; values = []; line_style = '-'; end
if nargin < 2; line_style = '-'; end
if isempty(values); values = []; end
if isempty(line_style); line_style = '-'; end

%% 1 - Plotting the horizontal lines at the y-values
axis_value = axis;
hold on;
for ii = 1:length(values)
    plot(axis_value(1:2), [values(ii),values(ii)], line_style);
end
hold off