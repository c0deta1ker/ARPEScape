function [xdat_crop, ydat_crop] = data_crop1D(xdat, ydat, xdat_lims, plot_result)
% [xdat_crop, ydat_crop] = data_crop1D(xdat, ydat, xdat_lims, plot_result)
%   This function crops 1D data along the x-axis (domain) and y-axis
%   (range). The user defines the limits along the x-axis. 
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   xdat:           [nX×1] array of the x-axis
%   -   ydat:           [nY×1] array of the y-axis
% 	-   xdat_lims:      [1×2] row vector of x-axis limits
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   xdat_crop:  	[mX×1] array of the cropped x-axis
%   -   ydat_crop:   	[mY×1] array of the cropped y-axis

%% Default parameters
if nargin < 3; xdat_lims=[min(xdat(:)), max(xdat(:))]; end
if nargin < 4; plot_result=0; end
if isempty(xdat_lims); xdat_lims=[min(xdat(:)), max(xdat(:))];  end
if isempty(plot_result); plot_result=0;  end
% - Sorting the cropping limits in ascending order
xdat_lims       = sort(xdat_lims);
%% 1 - Cropping the data
% - x-axis indices
[~, xIndxL]     = min(abs(xdat(:) - xdat_lims(1)));
[~, xIndxU]     = min(abs(xdat(:) - xdat_lims(2)));
x_indx          = [xIndxL xIndxU];
% - cropping x-axis
xdat_crop       = xdat(x_indx(1):x_indx(2));
ydat_crop       = ydat(x_indx(1):x_indx(2));
%% -- For Debugging
if plot_result == 1
    figure(); hold on;
    plot(xdat, ydat, 'k:', 'linewidth', 1);
    plot(xdat_crop, ydat_crop, 'k-', 'linewidth', 2);
    xline_n(xdat_lims, 'Color', [0 0 1], 'LineWidth', 1.5, 'Linestyle', '-');
    gca_props();
    xlabel(' X ', 'fontweight', 'bold');
    ylabel(' Y ', 'fontweight', 'bold');
    axis([min(xdat(:)), max(xdat(:)), min(ydat(:)), max(ydat(:))]);
end
end
