function [xdat_crop, ydat_crop, ddat_crop] = data_crop2D(xdat, ydat, ddat, xdat_lims, ydat_lims, plot_result)
% [xdat_crop, ydat_crop, ddat_crop] = data_crop2D(xdat, ydat, ddat, xdat_lims, ydat_lims, plot_result)
%   This function crops 2D data along the x- and y-axis. The user can
%   define the limits along the x- and y-axis to be cropped.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   xdat:           [1×nX] or [nY×nX] array of the x-axis data
%   -   ydat:           [nY×1] or [nY×nX] array of the y-axis data
%   -   ddat:           [nY×nX] array of the intensity data
% 	-   xdat_lims:      [1×2] row vector of x-axis limits
% 	-   ydat_lims:      [1×2] row vector of y-axis limits
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   xdat_crop:  	[mY×mX] array of the cropped x-axis
%   -   ydat_crop:   	[mY×mX] array of the cropped y-axis
%   -   ddat_crop:  	[mY×mX] array of the cropped intensity data

%% Default parameters
if nargin < 4; xdat_lims=[min(xdat(:)), max(xdat(:))]; ydat_lims=[min(ydat(:)), max(ydat(:))]; end
if nargin < 5; ydat_lims=[min(ydat(:)), max(ydat(:))]; end
if nargin < 6; plot_result=0; end
if isempty(xdat_lims); xdat_lims=[min(xdat(:)), max(xdat(:))];  end
if isempty(ydat_lims); ydat_lims=[min(ydat(:)), max(ydat(:))];  end
if isempty(plot_result); plot_result=0;  end
% - Sorting the cropping limits in ascending order
xdat_lims = sort(xdat_lims);
ydat_lims = sort(ydat_lims);
% - Initialising variables
xdat_crop = [];
ydat_crop = [];
ddat_crop = [];
%% 1 - Cropping the data
% - x-axis indices
if ~isempty(xdat)
    [~, xIndxL]     = min(abs(xdat(1,:) - xdat_lims(1)));
    [~, xIndxU]     = min(abs(xdat(1,:) - xdat_lims(2)));
    x_indx          = [xIndxL xIndxU];
end
% - y-axis indices
if ~isempty(ydat)
    [~, yIndxL]     = min(abs(ydat(:,1) - ydat_lims(1)));
    [~, yIndxU]     = min(abs(ydat(:,1) - ydat_lims(2)));
    y_indx          = [yIndxL yIndxU];
end
% - cropping x-axis
if ~isempty(xdat)
    if size(xdat, 1) > 1 && isempty(ydat); 	xdat_crop = xdat(:, x_indx(1):x_indx(2));
    elseif size(xdat, 1) > 1;               xdat_crop = xdat(y_indx(1):y_indx(2), x_indx(1):x_indx(2));
    else;                                   xdat_crop = xdat(1,x_indx(1):x_indx(2));
    end
end
% - cropping y-axis
if ~isempty(ydat)
    if size(ydat, 2) > 1 && isempty(xdat); 	ydat_crop = ydat(y_indx(1):y_indx(2), :);
    elseif size(ydat, 2) > 1;               ydat_crop = ydat(y_indx(1):y_indx(2), x_indx(1):x_indx(2));
    else;                                   ydat_crop = ydat(y_indx(1):y_indx(2),1);
    end
end
% - cropping the data
if isempty(ydat);       ddat_crop = ddat(:,x_indx(1):x_indx(2));
elseif isempty(xdat);   ddat_crop = ddat(y_indx(1):y_indx(2),:);
else;                   ddat_crop = ddat(y_indx(1):y_indx(2), x_indx(1):x_indx(2));
end
%% -- For Debugging
if plot_result == 1
    fig = figure(); 
    fig.Position(3) = 400*2; 
    fig.Position(4) = 400*0.9; 
    subplot(121); hold on;
    h1 = pcolor(xdat, ydat, ddat); set(h1,'EdgeColor','None','FaceColor','Flat');
    xpts = [xdat_lims(1), xdat_lims(2), xdat_lims(2), xdat_lims(1), xdat_lims(1)];
    ypts = [ydat_lims(2), ydat_lims(2), ydat_lims(1), ydat_lims(1), ydat_lims(2)];
    plot(xpts, ypts, 'Color', 'c', 'LineWidth', 1.5, 'Linestyle', '-');
    img_props();
    xlabel('$$ \bf  X $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Y $$', 'Interpreter', 'latex');
    axis([min(xdat(:)), max(xdat(:)), min(ydat(:)), max(ydat(:))]);
    subplot(122); hold on;
    h2 = pcolor(xdat_crop, ydat_crop, ddat_crop); set(h2,'EdgeColor','None','FaceColor','Flat');
    img_props();
    xlabel('$$ \bf  X $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  Y $$', 'Interpreter', 'latex');
    axis([min(xdat_crop(:)), max(xdat_crop(:)), min(ydat_crop(:)), max(ydat_crop(:))]);
end
end