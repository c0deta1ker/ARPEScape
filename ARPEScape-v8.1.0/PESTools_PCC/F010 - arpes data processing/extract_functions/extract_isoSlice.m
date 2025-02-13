function [sliceStr, fig] = extract_isoSlice(dataStr, isoslice_args, plot_results)
% [sliceStr, fig] = extract_isoSlice(dataStr, isoslice_args, plot_results)
%   This function plots the 3D ARPES data that is in the form 
%   D(X,Y,Z) into a single figure that includes the Eb(k) dispersion
%   and an IsoK or IsoE slice.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   [h=] ImData(X,Y,Z[,style]);
%   -   [DSlice,XSlice] = Slice(ACorr,ECorr,Data,xMode,Win) ;
%   -   [XR,YR,DataR] = Remap(XM,YM,Data);
%
%   IN:
%   -   dataStr:            data structure of the ARPES data.
%   -   isoslice_args:      1x4 cell of {scanIndex, isoType, isoWin, remap}.
%   -   plot_results:    	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -   sliceStr:    	new MATLAB data-structure.
%   -   fig:            figure output.

%% Default parameters
if nargin < 3; plot_results = 1; end
if nargin < 2; isoslice_args =  cell(1,4); plot_results = 1; end
if isempty(isoslice_args); isoslice_args = cell(1,4); end
if isempty(plot_results); plot_results = 1; end
% - Extracting the fields to be used with most recent processing
[xField, yField, zField, dField] = find_data_fields(dataStr);

%% - 1 - Initialising input parameters
scanIndex       = isoslice_args{1}; if isempty(scanIndex);  scanIndex   = 1; end            % Scan index to be cut
isoType         = isoslice_args{2}; if isempty(isoType);    isoType     = "IsoK"; end       % Type of cut to make. Either "IsoE" (horizontal) or "IsoK" (vertical)
isoWin          = isoslice_args{3}; if isempty(isoWin);     isoWin      = [-0.1, 0.1]; end  % Integration window of the cut to be made.
remap           = isoslice_args{4}; if isempty(remap);      remap       = 1; end            % Remap the slice onto a square grid? 1 or 0
% - Capping the window to the max/min of the ARPES data
if isempty(isoWin); Win = [-0.1, 0.1];
else; Win = sort(isoWin);
end
step_size = 0.025;
if isoType == "IsoE"
    if Win(1) < min(dataStr.(yField)(:)); Win(1) = min(dataStr.(yField)(:))+step_size; end
    if Win(2) > max(dataStr.(yField)(:)); Win(2) = max(dataStr.(yField)(:))-step_size; end
elseif isoType == "IsoK"
    if Win(1) < min(dataStr.(xField)(:)); Win(1) = min(dataStr.(xField)(:))+step_size; end
    if Win(2) > max(dataStr.(xField)(:)); Win(2) = max(dataStr.(xField)(:))-step_size; end
end
%% 2 - Extracting the Iso-Slice
if isoType == "IsoE"
    [DSlice, XSlice] = Slice(dataStr.(xField), dataStr.(yField), dataStr.(dField), 'IsoE', Win);
    % - Extracting the scan parameter variables
    if isfield(dataStr, 'kx')
        [~, indx] = min(abs(squeeze(dataStr.(yField)(:,1,1)) - mean(Win)));
        indx = ceil(indx); YSlice = squeeze(dataStr.(zField)(indx,:,:))';
    else; YSlice = squeeze(dataStr.(zField))'; YSlice = repmat(YSlice, [1, size(XSlice, 2)]);
    end
elseif isoType== "IsoK"
    [DSlice, YSlice] = Slice(dataStr.(xField), dataStr.(yField), dataStr.(dField), 'IsoK', Win);
    % - Extracting the scan parameter variables
    if isfield(dataStr, 'kx')
        [~, indx] = min(abs(squeeze(dataStr.(xField)(1,:,end)) - mean(Win)));
        indx = ceil(indx); XSlice = squeeze(dataStr.(zField)(:,indx,:));
    else; XSlice = squeeze(dataStr.(zField))'; XSlice = repmat(XSlice, [1, size(YSlice, 1)])';
    end
end
% - For a 3D Eb(kx,ky) scan type
if dataStr.Type == "Eb(kx,ky)";     scanAxis = dataStr.tltM;
% - For a 3D Eb(kx,kz) scan type
elseif dataStr.Type == "Eb(kx,kz)"; scanAxis = dataStr.hv;
% - For a 3D Eb(k,i) scan type
elseif dataStr.Type == "Eb(k,i)";   scanAxis = dataStr.index;
end
%% 3 - Re-mapping the Iso slices onto a square grid if required
if remap == 1
     [XSlice, YSlice, DSlice] = Remap(XSlice, YSlice, DSlice);
     YSlice = YSlice';
end
DSlice(isnan(DSlice)) = 0;
%% 4 - Appending data to MATLAB data structure
sliceStr                = struct();
sliceStr.isoslice_args  = isoslice_args;
sliceStr.ScanAxis       = scanAxis;
sliceStr.xField     	= string(xField);
sliceStr.FileName       = dataStr.FileName;
sliceStr.Type           = dataStr.Type;
sliceStr.IsoType        = isoType;
sliceStr.IsoWin         = Win;
sliceStr.IsoVal         = mean(Win);
sliceStr.IsoRan         = 0.5*range(Win);
sliceStr.IsoNum         = 1;
sliceStr.XSlice         = XSlice;
sliceStr.YSlice         = YSlice;
sliceStr.DSlice         = DSlice;
%% 5 - Plot the results
fig = [];
if plot_results == 1
    pp = plot_props();
    fig = figure(); 
    fig.Position(3) = 2*pp.fig5x4(1); 
    fig.Position(4) = 1*pp.fig5x4(2);
    if isfield(dataStr, 'kx');          ana_stage = "Algnd, Nrmlzed, Kconv";
    elseif isfield(dataStr, 'data');    ana_stage = "Algnd, Nrmlzed";
    elseif isfield(dataStr, 'eb');      ana_stage = "Algnd";
    else;                               ana_stage = "Raw";
    end
    set(fig, 'Name', sprintf("%s - %s: %s", dataStr.Type, isoType, ana_stage));
    %% 5.1 - Plotting the Eb(k) at the scan value
    subplot(1,2,1); hold on;
    % -- Plotting the ARPES Scan Image
    if isfield(dataStr, 'kx');                                  ImData(dataStr.(xField)(:,:,scanIndex), dataStr.(yField)(:,:,scanIndex), dataStr.(dField)(:,:,scanIndex));
    elseif isfield(dataStr, 'data') || isfield(dataStr, 'eb');  ImData(dataStr.(xField)(:,:,scanIndex), dataStr.(yField)(:,:,scanIndex), dataStr.(dField)(:,:,scanIndex));
    else;                                                       ImData(dataStr.(xField), dataStr.(yField), dataStr.(dField)(:,:,scanIndex));
    end
    % -- Plotting an outline over the integrated sliced region
    if isoType == "IsoE";       yline_n(Win, 'color', 'g', 'linewidth', 1);
    elseif isoType== "IsoK";    xline_n(Win, 'color', 'g', 'linewidth', 1);
    end
    % -- Formatting the figure
    img_props(); cbar_props([], 'Position', [0.06 0.10 0.015 0.175], 'YAxisLocation', 'left');
    axis([min(dataStr.(xField)(:)), max(dataStr.(xField)(:)),min(dataStr.(yField)(:)), max(dataStr.(yField)(:))]);
    title(sprintf(string(dataStr.FileName) + "; scan %i", scanIndex), 'interpreter', 'none');
    ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
    if isfield(dataStr, 'kx');      xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
    else;                           xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex');
    end
    %% 5.2 - Plotting the Iso-Slice
    subplot(1,2,2); hold on;
    % - Plotting the Iso slices
    ImData(XSlice, YSlice, DSlice);
    % -- Formatting the figure
    img_props(); cbar_props([], 'Position', [0.94 0.10 0.015 0.175], 'YAxisLocation', 'right');
    minC = min(DSlice(:)); maxC = max(DSlice(:)); clim([minC, maxC]);
    axis([min(XSlice(:)), max(XSlice(:)), min(YSlice(:)), max(YSlice(:))]);
    title(sprintf("%s; [%.2f,%.2f]", isoType, Win(1),Win(2)), 'interpreter', 'none');
    % - Re-labelling the axes depending on what slice is taken
    if isoType== "IsoE"
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
    elseif isoType== "IsoK"
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
            xlabel('$$ \bf  Index $$', 'Interpreter', 'latex');
        end
    end
    %% 5.3 - Extracting the curve formed by the Eb(k) scan slice being looked at
    if isoType== "IsoE"
        initIndx = floor(0.5*size(dataStr.(dField), 1));
        if isfield(dataStr, 'kx')
            x_curve = squeeze(dataStr.(xField)(initIndx,:,scanIndex));
            y_curve = squeeze(dataStr.(zField)(initIndx,:,scanIndex));
        elseif isfield(dataStr, 'data') || isfield(dataStr, 'eb')
            x_curve = squeeze(dataStr.(xField)(initIndx,:,scanIndex));
            y_curve = repmat(squeeze(dataStr.(zField)(scanIndex)), size(x_curve));
        else
            x_curve = squeeze(dataStr.(xField));
            y_curve = repmat(squeeze(dataStr.(zField)(scanIndex)), size(x_curve));
        end
        % - Fitting a 4th order polynomial to the slice curve
        x_interp = linspace(min(x_curve(:))-1e2, max(x_curve(:))+1e2, 1e3);
        y_interp = polyval(polyfit(x_curve,y_curve,4),x_interp);
    elseif isoType== "IsoK"
        initIndx = floor(0.5*size(dataStr.(dField), 2));
        if isfield(dataStr, 'kx')
            y_curve = squeeze(dataStr.(yField)(:,initIndx,scanIndex));
            x_curve = squeeze(dataStr.(zField)(:,initIndx,scanIndex));
            x_interp = linspace(min(x_curve(:))-5, max(x_curve(:))+5, 1e3);
            y_interp = polyval(polyfit(x_curve,y_curve,1),x_interp);
        elseif isfield(dataStr, 'data') || isfield(dataStr, 'eb')
            y_curve = squeeze(dataStr.(yField)(:,initIndx,scanIndex));
            x_curve = repmat(squeeze(dataStr.(zField)(scanIndex)), size(y_curve));
            x_interp = ones(size(x_curve))*mean(x_curve(:));
            y_interp = linspace(-1e3, 1e3, size(x_interp, 1));
        else
            y_curve = squeeze(dataStr.(yField));
            x_curve = repmat(squeeze(dataStr.(zField)(scanIndex)), size(y_curve));
            x_interp = ones(size(x_curve))*mean(x_curve(:));
            y_interp = linspace(-1e3, 1e3, size(x_interp, 1));
        end
    end
    % Plotting the scan slice curve
    plot(x_interp, y_interp, 'Color', [0 0 1], 'LineWidth', 1., 'Linestyle', '-');
end
end