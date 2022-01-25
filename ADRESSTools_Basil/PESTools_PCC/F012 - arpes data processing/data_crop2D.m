function [xDat_crop, yDat_crop, dDat_crop] = data_crop2D(xDat, yDat, dDat, xDat_lims, yDat_lims)
% [xDat_crop, yDat_crop, dDat_crop] = data_crop2D(xDat, yDat, dDat, xDat_lims, yDat_lims)
%   This function crops the ARPES x-, y- and z-independent variables over
%   a given range. The function crops the most recently processed data and
%   ensures consistency across both the independent variables and data
%   matrix.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%
%   IN:
%   -   xDat:           [N x M] array of the x-axis (theta/kx)
%   -   yDat:           [N x M] array of the y-axis (eb)
%   -   dDat:           [N x M] array of the ARPES data (intensity)
% 	-   xDat_lims:      [1 x 2] row vector of x-axis limits
% 	-   yDat_lims:      [1 x 2] row vector of y-axis limits
%
%   OUT:
%   -   xDat_crop:  	[A x B] array of the cropped x-axis (theta/kx)
%   -   yDat_crop:   	[A x B] array of the cropped y-axis (eb)
%   -   dDat_crop:  	[A x B] array of the cropped ARPES data (intensity)

%% Default parameters
maxLim = 1e4;
if nargin < 4; xDat_lims=[-1 1]*maxLim; yDat_lims=[-1 1]*maxLim; end
if nargin < 5; yDat_lims=[-1 1]*maxLim; end
if isempty(xDat_lims); xDat_lims=[-1 1]*maxLim;  end
if isempty(yDat_lims); yDat_lims=[-1 1]*maxLim;  end
% - Sorting the cropping limits in ascending order
xDat_lims = sort(xDat_lims);
yDat_lims = sort(yDat_lims);
% - Initialising variables
xDat_crop = [];
yDat_crop = [];
dDat_crop = [];

%% 1 - Cropping the data
% - x-axis indices
if ~isempty(xDat)
    [~, xIndxL]     = min(abs(xDat(1,:) - xDat_lims(1)));
    [~, xIndxU]     = min(abs(xDat(1,:) - xDat_lims(2)));
    x_indx          = [xIndxL xIndxU];
end

% - y-axis indices
if ~isempty(yDat)
    [~, yIndxL]     = min(abs(yDat(:,1) - yDat_lims(1)));
    [~, yIndxU]     = min(abs(yDat(:,1) - yDat_lims(2)));
    y_indx          = [yIndxL yIndxU];
end

% - cropping x-axis
if ~isempty(xDat)
    if size(xDat, 1) > 1 && isempty(yDat); 	xDat_crop = xDat(:, x_indx(1):x_indx(2));
    elseif size(xDat, 1) > 1;               xDat_crop = xDat(y_indx(1):y_indx(2), x_indx(1):x_indx(2));
    else;                                   xDat_crop = xDat(1,x_indx(1):x_indx(2));
    end
end
% - cropping y-axis
if ~isempty(yDat)
    if size(yDat, 2) > 1 && isempty(xDat); 	yDat_crop = yDat(y_indx(1):y_indx(2), :);
    elseif size(yDat, 2) > 1;               yDat_crop = yDat(y_indx(1):y_indx(2), x_indx(1):x_indx(2));
    else;                                   yDat_crop = yDat(y_indx(1):y_indx(2),1);
    end
end
% - cropping the data
if isempty(yDat)
    dDat_crop = dDat(:,x_indx(1):x_indx(2));
elseif isempty(xDat)
    dDat_crop = dDat(y_indx(1):y_indx(2),:);
else
    dDat_crop = dDat(y_indx(1):y_indx(2), x_indx(1):x_indx(2));
end
end
