function [cutStr, fig] = extract_isoCut(dataStr, isocut_args, plot_results)
% [cutStr, fig] = extract_isoCut(dataStr, isocut_args, plot_results)
%   This function plots a EDC/MDC cut through the ARPES data in 
%   D(X,Y,cutIndx) into a single figure that includes the Eb(k) 
%   dispersion.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   [h=] ImData(X,Y,Z[,style])
%   -   [XCut,DCut] = Cut(ACorr,ECorr,Data,xMode,Win)
%
%   IN:
%   -   dataStr:      	data structure of the ARPES data.
%   -   isocut_args:  	1x3 cell of {scanIndex, isoType, isoWin}.
%   -   plot_results:  	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -   cutStr:         new MATLAB data-structure.
%   -   fig:            figure output.

%% Default parameters
if nargin < 3; plot_results = 1; end
if nargin < 2; isocut_args = cell(1,3); plot_results = 1; end 
if isempty(isocut_args); isocut_args = cell(1,3); end 
if isempty(plot_results); plot_results = 1; end
pp = plot_props();
fig = [];
              
%% - 1 - Initialising input parameters
scanIndex       = isocut_args{1}; if isempty(scanIndex); scanIndex = 1; end         % Scan index to be cut
isoType         = isocut_args{2}; if isempty(isoType); isoType = "EDC"; end         % Type of cut to make. Either "MDC" (horizontal) or "EDC" (vertical)
isoWin          = isocut_args{3}; if isempty(isoWin); isoWin = [-0.1, 0.1]; end     % Integration window of the cut to be made.
% - Capping the window to the max/min of the ARPES data
if isempty(isoWin); Win = [-0.1, 0.1];
else; Win = sort(isoWin);
end
% - Extracting the x-axis, y-axis and data
if isfield(dataStr, 'IsoType')
    scanIndex   = 1;
    xField      = dataStr.xField;
    if dataStr.IsoType == "IsoE" || dataStr.IsoType == "IsoK"
        XDat         = dataStr.XSlice;
        YDat         = dataStr.YSlice;
        DDat         = dataStr.DSlice;
    elseif dataStr.IsoType == "Scan"
        XDat         = dataStr.XScan;
        YDat         = dataStr.YScan;
        DDat         = dataStr.DScan;
    end
else
    [xField, yField, ~, dField] = find_data_fields(dataStr);
    if isfield(dataStr, 'kx') || isfield(dataStr, 'data') || isfield(dataStr, 'eb')
        XDat         = dataStr.(xField)(:,:,scanIndex);
        YDat         = dataStr.(yField)(:,:,scanIndex);
        DDat         = dataStr.(dField)(:,:,scanIndex);
    else
        XDat         = dataStr.(xField);
        YDat         = dataStr.(yField);
        DDat         = dataStr.(dField)(:,:,scanIndex);
    end
end
% - Capping the window limits
step_size = 0.025;
if isoType == "MDC"
    if Win(1) < min(YDat(:)); Win(1) = min(YDat(:)) + step_size; end
    if Win(2) > max(YDat(:)); Win(2) = max(YDat(:)) - step_size; end
elseif isoType == "EDC"
    if Win(1) < min(XDat(:)); Win(1) = min(XDat(:)) + step_size; end
    if Win(2) > max(XDat(:)); Win(2) = max(XDat(:)) - step_size; end
end
%% 2 - Extracting the Iso cuts
if isoType == "MDC"
    [XCut, DCut] = Cut(XDat, YDat, DDat, 'mdc', Win);
elseif isoType== "EDC"
    [XCut, DCut] = Cut(XDat, YDat, DDat, 'edc', Win);
end

%% 3.0 - Plotting the figure if required
if plot_results == 1
    fig = figure(); 
    fig.Position(1) = pp.figpos(1);
    fig.Position(2) = pp.figpos(2);
    if isoType == "MDC"
        fig.Position(3) = pp.fig5x4(1); 
        fig.Position(4) = 1.5*pp.fig5x4(2);
        subplot(3,1,2:3); hold on;
    elseif isoType== "EDC"
        fig.Position(3) = 1.5*pp.fig5x4(1); 
        fig.Position(4) = pp.fig5x4(2);
        subplot(1,3,1:2); hold on;
    end
    % - For an Eb(k) scan type
    if dataStr.Type == "Eb(k)"
        % -- Plotting the Eb(k) figure object
        if isfield(dataStr, 'kx'); set(fig, 'Name', 'Iso-Cut - Eb(k): A, N, K');
        elseif isfield(dataStr, 'data'); set(fig, 'Name', 'Iso-Cut - Eb(k): A, N');
        elseif isfield(dataStr, 'eb'); set(fig, 'Name', 'Iso-Cut - Eb(k): A');
        else; set(fig, 'Name', 'Iso-Cut - Eb(k): Raw');
        end
    % - For a 3D Eb(kx,ky) scan type
    elseif dataStr.Type == "Eb(kx,ky)"
        % -- Plotting the Eb(kx,ky) figure object
        if isfield(dataStr, 'kx'); set(fig, 'Name', 'Iso-Cut - Eb(kx,ky): A, N, K');
        elseif isfield(dataStr, 'data'); set(fig, 'Name', 'Iso-Cut - Eb(kx,ky): A, N');
        elseif isfield(dataStr, 'eb'); set(fig, 'Name', 'Iso-Cut - Eb(kx,ky): A');
        else; set(fig, 'Name', 'Iso-Cut - Eb(kx,ky): Raw');
        end
    % - For a 3D Eb(kx,kz) scan type
    elseif dataStr.Type == "Eb(kx,kz)"
        % -- Plotting the Eb(kx,kz) figure object
        if isfield(dataStr, 'kx'); set(fig, 'Name', 'Iso-Cut - Eb(kx,kz): A, N, K');
        elseif isfield(dataStr, 'data'); set(fig, 'Name', 'Iso-Cut - Eb(kx,kz): A, N');
        elseif isfield(dataStr, 'eb'); set(fig, 'Name', 'Iso-Cut - Eb(kx,kz): A');
        else; set(fig, 'Name', 'Iso-Cut - Eb(kx,kz): Raw');
        end
    % - For a 3D Eb(k,i) scan type
    elseif dataStr.Type == "Eb(k,i)"
        % -- Plotting the Eb(kx,kz) figure object
        if isfield(dataStr, 'kx'); set(fig, 'Name', 'Iso-Cut - Eb(k,i): A, N, K');
        elseif isfield(dataStr, 'data'); set(fig, 'Name', 'Iso-Cut - Eb(k,i): A, N');
        elseif isfield(dataStr, 'eb'); set(fig, 'Name', 'Iso-Cut - Eb(k,i): A');
        else; set(fig, 'Name', 'Iso-Cut - Eb(k,i): Raw');
        end
    end
    %% 3.1 - Plotting the Eb(k) at the scan value
    % -- The plot depends on how far in the analysis 
    ImData(XDat, YDat, DDat);
    img_props([], string(xField));
    axis([min(XDat(:)), max(XDat(:)),min(YDat(:)), max(YDat(:))]);
    % -- Plotting an outline over the integrated sliced region
    if isoType == "MDC"
        colbar_props();
        patch([-1e3, 1e3, 1e3, -1e3, -1e3], [Win(1), Win(1), Win(2), Win(2), Win(1)],...
            [0 1 0], 'linewidth', 1, 'facealpha', 0, 'edgecolor', [0 1 0]);
    elseif isoType== "EDC"
        colbar_props([0.05 0.125 0.025 0.175], 'left');
        patch([Win(1), Win(1), Win(2), Win(2), Win(1)],[-1e3, 1e3, 1e3, -1e3, -1e3],...
            [0 1 0], 'linewidth', 1, 'facealpha', 0, 'edgecolor', [0 0 1]);
    end
    % -- Adding title to the figure
    title(sprintf(string(dataStr.FileName) + "; scan %i", scanIndex), 'interpreter', 'none', 'fontsize', 9);
% (A) -- If the data to be saved is an IsoSlice / IsoCut /IsoScan form
    if isfield(dataStr, 'IsoType')
        if dataStr.IsoType== "IsoE"
            if dataStr.Type == "Eb(kx,ky)"
                if isfield(dataStr, 'kx'); xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  \tau (^{\circ}) $$', 'Interpreter', 'latex');
                end
            elseif dataStr.Type == "Eb(kx,kz)"
                if isfield(dataStr, 'kx'); xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
                end
            elseif dataStr.Type == "Eb(k,i)"
                if isfield(dataStr, 'kx'); xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  Index $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  Index $$', 'Interpreter', 'latex'); ylabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
                end
            end
        elseif dataStr.IsoType== "IsoK"
            ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex'); 
            if dataStr.Type == "Eb(kx,ky)"
                if isfield(dataStr, 'kx'); xlabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex');
                else; xlabel('$$ \bf  \tau (^{\circ}) $$', 'Interpreter', 'latex');
                end
            elseif dataStr.Type == "Eb(kx,kz)"
                if isfield(dataStr, 'kx'); xlabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
                end
            elseif dataStr.Type == "Eb(k,i)"
                if isfield(dataStr, 'kx'); xlabel('$$ \bf  Index $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  Index $$', 'Interpreter', 'latex');
                end
            end
% (B) -- Else, for the usual ARPES data structure
        else
            img_props([], string(dataStr.xField));
        end
    end
    %% 3.2 - Plotting the Iso-Cut
    if isoType == "MDC"
        subplot(3,1,1); hold on;
        plot(XCut, DCut, 'k.-', 'color', [0 1 0], 'LineWidth', 1.5);
        gca_props();
        ax = gca;
        ax.XAxisLocation = 'top';               % 'bottom' | 'top' | 'origin'
        ax.YAxisLocation = 'left';              % 'left' | 'right' | 'origin'
        % Axis labels and limits
        if xField == "raw_tht" || xField == "tht"
            xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex');
        elseif xField == "kx"
            xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
        end
        ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
        % - Plotting the x- and y-axes
        xl = [min(XDat(:)), max(XDat(:))]; yl = ylim;
        axis([xl(1), xl(2), yl(1), 1.05*yl(2)]);
    elseif isoType == "EDC"
        subplot(1,3,3); hold on;
        plot(DCut, XCut, 'k.-', 'color', [0 0 1], 'LineWidth', 1.5);
        gca_props();
        ax = gca;
        ax.XAxisLocation = 'bottom';            % 'bottom' | 'top' | 'origin'
        ax.YAxisLocation = 'right';             % 'left' | 'right' | 'origin'
        % - Axis labels and limits
        ylabel('$$ \bf  E_{B} (eV) $$', 'Interpreter', 'latex');
        xlabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
        % - Plotting the x- and y-axes
        xl = xlim; yl = [min(YDat(:)), max(YDat(:))];
        axis([xl(1), 1.05*xl(2), yl(1), yl(2)]);
    end
end

%% 4 - Appending data to MATLAB data structure
cutStr              = struct();
cutStr.isocut_args  = isocut_args;
cutStr.xField     	= string(xField);
cutStr.FileName     = dataStr.FileName;
cutStr.Type         = dataStr.Type;
cutStr.IsoType   	= isoType;
cutStr.IsoVal   	= mean(Win);
cutStr.IsoWin   	= 0.5*range(Win);
cutStr.IsoNum       = 1;
cutStr.XCut         = XCut;
cutStr.DCut         = DCut;
end