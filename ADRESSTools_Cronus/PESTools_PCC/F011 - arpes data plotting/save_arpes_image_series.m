function save_arpes_image_series(dataStr, viewVideo_args)
% save_arpes_image_series(dataStr, viewVideo_args)
%   This function plots the Eb(kx) ARPES image as a function
%   of the scan parameter, with the addition of a zoomed-in area 
%   of interest as an image series. You can also use this function to save
%   all the ARPES scan images as an image series, allowing you to browse
%   through all the ARPES scan images easily.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   gca_properties(type)
%   -   [h=] ImData(X,Y,Z[,style])
%
%   IN:
%   -   dataStr:        data structure of the ARPES data.
%   -   viewVideo_args: 1x7 cell of {zLims, xLims, yLims, interpolate, ROI, xLimsROI, yLimsROI}.
%
%   OUT:
%   -   figure output.

%% Default parameters
pp = plot_props();
if nargin < 2; viewVideo_args = cell(1,7); end
if isempty(viewVideo_args); viewVideo_args = cell(1,7); end
% - Extracting the fields to be used with most recent processing
[xField, yField, zField, dField] = find_data_fields(dataStr);

%% - 1 - Initialising the defined parameters
zLims           = viewVideo_args{1}; if isempty(zLims); zLims = []; end
xLims           = viewVideo_args{2}; if isempty(xLims); xLims = []; end
yLims           = viewVideo_args{3}; if isempty(yLims); yLims = []; end
interpolate     = viewVideo_args{4}; if isempty(interpolate); interpolate = 0; end
ROI             = viewVideo_args{5}; if isempty(ROI); ROI = 0; end
xLimsROI        = viewVideo_args{6}; if isempty(xLimsROI); xLimsROI = []; end
yLimsROI        = viewVideo_args{7}; if isempty(yLimsROI); yLimsROI = []; end

%% - 2 - Initialising the figure
% - Choosing where to save the data
filter = {'*.*'};
[save_filename, save_filepath] = uiputfile(filter, 'Save the image series...', dataStr.FileName);
fName = string(string(save_filepath) + string(save_filename));
% If Cancel is pressed, then return nothing
if isequal(save_filepath,0) || isequal(save_filename,0); return; end
% - Initialising the data structures
dataStr = data_crop3D(dataStr, [], [], zLims);
% - Initialising the frames
nFrames = size(dataStr.(dField), 3);

%% - 3 - Running the Eb(k) vs scan parameter frame-by-frame
for f = 1:nFrames
    fig = figure(); 
    fig.Position(3) = 1.25*pp.fig5x4(1); 
    fig.Position(4) = 1.25*pp.fig5x4(2);
    % - Initialising variables
    ifName = fName + "_" + string(f);
    % - Extracting the scan parameter for the given frame
    if dataStr.Type == "Eb(kx,ky)" 
        if isfield(dataStr, 'kx'); scanVal = sprintf('$$ \\bf ky_{%i} = %.3f \\AA^{-1}$$', f, dataStr.(zField)(1,1,f));
        else; scanVal = sprintf('$$ \\bf \\tau_{%i} = %.2f^{\\circ}$$', f, dataStr.(zField)(f));
        end
    elseif dataStr.Type == "Eb(kx,kz)" 
        if isfield(dataStr, 'kx'); scanVal = sprintf('$$ \\bf kz_{%i} = %.3f \\AA^{-1}$$', f, dataStr.(zField)(1,1,f));
        else; scanVal = sprintf('$$ \\bf hv_{%i} = %.2f eV $$', f, dataStr.(zField)(f));
        end
    end
    % - Extracting the full data
    if isfield(dataStr, 'kx')
        xDat = dataStr.(xField)(:,:,f);
        yDat = dataStr.(yField)(:,:,f);
        dDat = dataStr.(dField)(:,:,f);
    elseif isfield(dataStr, 'data') || isfield(dataStr, 'eb')
        xDat = dataStr.(xField)(:,:,f);
        yDat = dataStr.(yField)(:,:,f);
        dDat = dataStr.(dField)(:,:,f);
    else
        xDat = dataStr.(xField);
        yDat = dataStr.(yField);
        dDat = dataStr.(dField)(:,:,f);
    end
    %% 2.1 - Plotting the full Eb(k) image first
    hold on;
    if interpolate == 1;        ImData(xDat, yDat, dDat, 'interp');
    else;                       ImData(xDat, yDat, dDat);
    end
    % - 2.2 - Formatting the figure
    minC = min(min(min(dataStr.(dField)(:,:,f))));
    maxC = max(max(max(dataStr.(dField)(:,:,f))));
    caxis([minC, maxC]);
    axis([xLims, yLims]);
    img_props([], string(xField));
    colbar_props();
    title(dataStr.FileName, 'interpreter', 'none', 'fontsize', 12);
    % - 2.3 - Adding text for the scan parameter
    annotation('textbox', [0.22 0.2 0.17 0.06], 'String',scanVal, 'FitBoxToText','on',...
        'color', [0 0 0], 'fontsize', 13, 'backgroundcolor', [1 1 1], 'facealpha', 0.80,...
        'linewidth', 2, 'horizontalalignment', 'center', 'verticalalignment', 'middle',...
        'interpreter', 'latex');
    %% 3 - Plot the ROI if required
    if ROI == 1
        %% 3.1 - Plotting the region of interest that is zoomed in on
        x = [xLimsROI(1), xLimsROI(2), xLimsROI(2), xLimsROI(1), xLimsROI(1)];
        y = [yLimsROI(1), yLimsROI(1), yLimsROI(2), yLimsROI(2), yLimsROI(1)] ;
        plot(x, y, 'LineWidth', 1.5, 'Color', 'y', 'linestyle', '-');
        %% 3.2 - Plotting the cropped Eb(k) on top
        % - Creating a new inset axes to the figure
        if mean(xLimsROI) < 0; new_ax = axes('Position',[.67 .655 .25 .25]);
        else; new_ax = axes('Position',[.12 .655 .25 .25]);
        end
        % Cropping the data
        [xDat_roi, yDat_roi, dDat_roi] = data_crop2D(xDat, yDat, dDat, xLimsROI, yLimsROI);
        % - Plotting the cropped ARPES data
        if interpolate == 1;        ImData(xDat_roi, yDat_roi, dDat_roi, 'interp');
        else;                       ImData(xDat_roi, yDat_roi, dDat_roi);
        end
        % - 3.3 - Formatting the figure
        minC = min(dDat_roi(:));
        maxC = max(dDat_roi(:));
        caxis([minC, maxC]);
        axis([min(xDat_roi(:)), max(xDat_roi(:)), min(yDat_roi(:)), max(yDat_roi(:))]);
        img_props([], string(xField));
        xlabel(''); ylabel('');
        % Other properties
        ax = gca;  
        ax.FontSize = 10; 
        ax.TickDir = 'both';
        ax.Color = [1 1 1]; 
        ax.LineWidth = 0.75;
        ax.Box = 'on'; 
        ax.Layer = 'Top';
        new_ax.XColor = [1 1 1]; 
        new_ax.YColor = [1 1 1]; 
        pbaspect([1 1 1]);
        if mean(xLimsROI) < 0
            new_ax.XAxisLocation = 'bottom';
            new_ax.YAxisLocation = 'left';
        else
            new_ax.XAxisLocation = 'bottom';
            new_ax.YAxisLocation = 'right';
        end
    end
    %% SAVE THE FIGURE
    print(char(ifName), '-dpng', '-r500');
end
end