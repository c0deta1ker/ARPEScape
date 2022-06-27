function fig = view_arpes_isoSlice(sliceStr, plot_args, bzOverlay, bz_args)
% fig = view_arpes_isoSlice(sliceStr, plot_args, bzOverlay, bz_args)
%   This function can be used to plot a very nice, publication quality 
%   IsoE image, along with the BZ overlay. Use this as a general function
%   to plot good IsoE images, where many plot arguments can be controlled.
%   The input data should be a sliceStr, that comes directly from the 
%   'extract_isoSlice()' function.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   sliceStr:     	data structure of the sliced Iso-E ARPES data.
%   -   plot_args:     	1x10 cell of {cMap, xLims, yLims, cLims, axCol, axisOn, axisEqual, gridOn, bzOn, interpOn}; 
%   -   bzOverlay:    	data structure of the sliced Iso-E ARPES data.
%   -   bz_args:        1x3 cell of {bz_offset [xShift,yShift], bz_col [1,1,1], bz_lwidth [1.5]}
%
%   OUT: 
%   -   fig:            MATLAB figure object with the ARPES data plotted.

%% Default parameters
if nargin < 4; bz_args = {[0, 0], [1 1 1], 1.5};  end
if nargin < 3; bzOverlay = []; bz_args = {[0, 0], [1 1 1], 1.5};  end
if nargin < 2; plot_args = cell(1,10); bzOverlay = []; bz_args = {[0, 0], [1 1 1], 1.5};  end
if isempty(bz_args);    bz_args = {[0, 0], [1 1 1], 1.5}; end
if isempty(bzOverlay);  bzOverlay = []; end
if isempty(plot_args);  plot_args = cell(1,10); end

%% - 1 - Initialising the parameters
% - Extracting the plot arguments
cMap        = plot_args{1};     if isempty(cMap);	cMap  	= "hot"; end    % colormap
xLims       = plot_args{2};     if isempty(xLims);  xLims   = [min(sliceStr.XSlice(:)), max(sliceStr.XSlice(:))]; end	% x-axis limits
yLims       = plot_args{3};     if isempty(yLims);  yLims   = [min(sliceStr.YSlice(:)), max(sliceStr.YSlice(:))]; end	% y-axis limits
cLims       = plot_args{4};     if isempty(cLims);  cLims   = [0, 1]; end           % color limits
axCol       = plot_args{5};     if isempty(axCol);  axCol   = [0, 0, 0]; end        % x- and y-axis colors
axisOn      = plot_args{6};     if isempty(axisOn); axisOn  = 1; end                % show axis lines?
axisEqual   = plot_args{7};     if isempty(axisEqual);	axisEqual   = 0; end        % plot the axis as equals
gridOn      = plot_args{8};     if isempty(gridOn);     gridOn      = 1; end        % show Grid?
bzOn        = plot_args{9};     if isempty(bzOn);    	bzOn        = 0; end        % show BZ?
interpOn 	= plot_args{10};    if isempty(interpOn);   interpOn    = 0; end        % interpolate the data?
% - Initialising the plot / figure properties
pp  = plot_props();
fig = figure(); hold on;
fig.Position(3) = pp.fig5x4(1);
fig.Position(4) = pp.fig5x4(2);

%% - 2 - Plotting the iso-slices
% 2.1 - Plotting the Iso-Slice plot
if interpOn == 1;   ImData(sliceStr.XSlice, sliceStr.YSlice, sliceStr.DSlice, 'interp');
else;               ImData(sliceStr.XSlice, sliceStr.YSlice, sliceStr.DSlice);
end

% - Initial formatting of the axis
ax = gca;
% -- Font properties
ax.FontName         = 'Helvetica'; 
ax.FontWeight       = 'normal'; 
ax.FontSize         = 12;
% -- Tick properties
ax.TickLabelInterpreter = 'latex';
ax.XMinorTick       = 'off'; 
ax.YMinorTick       = 'off';
ax.TickDir          = 'in';
ax.TickLength       = [0.02 0.025];
ax.XColor           = axCol; 
ax.YColor           = axCol;
% -- Ruler properties
ax.XAxisLocation    = 'bottom';            % 'bottom' | 'top' | 'origin'
ax.YAxisLocation    = 'left';              % 'left' | 'right' | 'origin'
% -- Box Styling properties
ax.Color            = [0, 0, 0];
ax.LineWidth        = 2.0;
ax.Box              = 'on'; 
ax.Layer            = 'Top';
% -- Adding title
title(sprintf(string(sliceStr.FileName) +...
    "; %s; [%.2f,%.2f]", sliceStr.IsoType, sliceStr.IsoWin(1), sliceStr.IsoWin(2)), 'interpreter', 'none', 'fontsize', 8);

% - Initial formatting of the figure
fig.Color           = [1 1 1]; 
fig.InvertHardcopy  = 'off';

% 2.2 - Plotting the Brilluoin zone overlayer
if bzOn == 1 && ~isempty(bzOverlay)
    bzone_props(bzOverlay, bz_args);
end

% 2.3  - Applying color and axes limits, grid- and axes-lines
if axisEqual == 1; axis equal; end
axis([xLims, yLims]); 
maxC = max(sliceStr.DSlice(:));
caxis([cLims(1), cLims(2)].*maxC);
if axisOn == 1
    line([0 0], [-1e5, 1e5], 'Color', axCol, 'LineWidth', 1, 'Linestyle', '--');
    line([-1e5, 1e5], [0 0], 'Color', axCol, 'LineWidth', 1, 'Linestyle', '--');
end
if gridOn == 1; grid on; end

% - Formatting the figure
colormap(cMap);
colbar_props();
% - Re-labelling the axes depending on what slice is taken
if sliceStr.IsoType== "IsoE"
    if sliceStr.Type == "Eb(kx,ky)"
        if sliceStr.xField == "kx"; xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex'); 
        else; xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  \tau (^{\circ}) $$', 'Interpreter', 'latex');
        end
    elseif sliceStr.Type == "Eb(kx,kz)"
        if sliceStr.xField == "kx"; xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex'); 
        else; xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
        end
    end
elseif sliceStr.IsoType== "IsoK"
    ylabel('$$ \bf  E_B (eV) $$', 'Interpreter', 'latex'); 
    if sliceStr.Type == "Eb(kx,ky)"
        if sliceStr.xField == "kx"; xlabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex');
        else; xlabel('$$ \bf  \tau (^{\circ}) $$', 'Interpreter', 'latex');
        end
    elseif sliceStr.Type == "Eb(kx,kz)"
        if sliceStr.xField == "kx"; xlabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex'); 
        else; xlabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
        end
    end
end

end