function fig = view_arpes_data_3D(dataStr, slice_args, varargin)
% fig = view_arpes_data_3D(dataStr, slice_args, varargin)
%   Use this as the general function to view all 3D ARPES data, both during
%   and after the processing of the data. The 3D figure is plotted with 3
%   panels; (1) IsoScan, (2) IsoE and (3) IsoK slices through the data along 
%   each one of the data dimensions. This function ensures that the most recent 
%   data processing is applied.
%
%   IN:
%   -   dataStr:    data structure of the ARPES data.
%   -   slice_args: 1x3 cell of slice values for {IsoScan, IsoE, IsoK}.
%   -   varargin:   AxProps arguments: Axes properties control the appearance and behavior of an Axes object, specified as a string scalar or character vector.
%
%   OUT:
%   -   fig:        Final MATLAB figure object. 

%% Default parameters
if nargin < 2; slice_args = cell(1,3); end
if nargin < 3; varargin = {}; end
if isempty(slice_args); slice_args = cell(1,3); end
if isempty(varargin); varargin = {}; end
% - Extracting the fields to be used with most recent processing
[xField, yField, zField, dField] = find_data_fields(dataStr);
pp = plot_props();
%% - 1 - Validity check on the input parameters
% -- Checking the length of input arguments
if length(slice_args) == 1; slice_args = {slice_args{1}, [], []};
elseif length(slice_args) == 2; slice_args = {slice_args{1}, slice_args{2}, []};
elseif length(slice_args) > 3; error('Number of slice arguments exceeds number of varables. Index must not exceed 3.');
end
% -- Loading in initial values
isoScan_val = slice_args{1}; if isempty(isoScan_val) || isnan(isoScan_val); isoScan_val = mean(dataStr.(zField)(:)); end     % Scan value to be sliced
isoE_val    = slice_args{2}; if isempty(isoE_val) || isnan(isoE_val);       isoE_val    = mean(dataStr.(yField)(:)); end     % IsoE value to be sliced
isoK_val    = slice_args{3}; if isempty(isoK_val) || isnan(isoK_val);       isoK_val    = mean(dataStr.(xField)(:)); end     % IsoK value to be sliced
% -- Padding values to the min/max values of the ARPES data
if isoK_val < min(dataStr.(xField)(:)); isoK_val = min(dataStr.(xField)(:)); end
if isoK_val > max(dataStr.(xField)(:)); isoK_val = max(dataStr.(xField)(:)); end
if isoE_val < min(dataStr.(yField)(:)); isoE_val = min(dataStr.(yField)(:)); end
if isoE_val > max(dataStr.(yField)(:)); isoE_val = max(dataStr.(yField)(:)); end
if isoScan_val < min(dataStr.(zField)(:)); isoScan_val = min(dataStr.(zField)(:)); end
if isoScan_val > max(dataStr.(zField)(:)); isoScan_val = max(dataStr.(zField)(:)); end
%% - 2 - Defining the relevant data for plotting
% -- Extracting the variables from the arpes data structure
XX  = dataStr.(xField);
YY  = flipud(dataStr.(yField));
ZZ  = dataStr.(zField);
DD  = dataStr.(dField); DD = flipud(DD);
% -- Squeezing the data into a matrix
D1  = double(squeeze(DD)); 
DIM = size(D1); 
% -- Extracting the bisecting planes for 3D data
% --- First determine the zField index
if string(zField) == "hv" || string(zField) == "tltM" || string(zField) == "index"; [~, indxZ] = min(abs(ZZ(1,:) - isoScan_val));
else;                                                                               [~, indxZ] = min(abs(ZZ(round(DIM(1)/2),round(DIM(2)/2),:) - isoScan_val));
end
% --- Determine the xField index
if string(xField) == "raw_tht";     [~, indxX] = min(abs(XX(1,:) - isoK_val));
else;                               [~, indxX] = min(abs(XX(round(DIM(2)/2),:,indxZ) - isoK_val));
end
% --- Determine the yField index
if string(yField) == "raw_eb";     [~, indxY] = min(abs(YY(:,1) - isoE_val));
else;                              [~, indxY] = min(abs(YY(:,indxX,indxZ) - isoE_val));
end
%% - 3 - Plotting the final 2D Figure
fig = figure(); 
fig.Position(3) = 2*pp.fig5x4(1); 
fig.Position(4) = pp.fig5x4(2);
% -- Set figure name
if isfield(dataStr, 'kx');          ana_stage = "Algnd, Nrmlzed, Kconv";
elseif isfield(dataStr, 'data');    ana_stage = "Algnd, Nrmlzed";
elseif isfield(dataStr, 'eb');      ana_stage = "Algnd";
else;                               ana_stage = "Raw";
end
set(fig, 'Name', sprintf("%s - %s", dataStr.Type, ana_stage));
% -- Creating a tiled axis
t = tiledlayout(1,3);
t.TileSpacing = 'compact';
t.Padding = 'compact';
title(t,dataStr.FileName,'FontWeight','bold','interpreter', 'none');
%% -- Plotting the IsoScan ARPES Image
nexttile;
if ndims(YY) == 2;  imagesc('XData', XX, 'YData', YY, 'CData', D1(:,:,indxZ), [min(D1(:)) max(D1(:))]);
else;               surf(squeeze(XX(:,:,indxZ)),squeeze(YY(:,:,indxZ)), squeeze(D1(:,:,indxZ))); view(2); shading flat; clim([min(D1(:)) max(D1(:))]);
end
img_props(0, 'Color', [1 1 1], varargin{:}); colormap(hot); 
title(sprintf('coronal (XY): IsoScan = %.2f', isoScan_val)); 
axis([min(XX(:)), max(XX(:)), min(YY(:)), max(YY(:))]);
xline(isoK_val, 'color', 'b', 'LineWidth', 1.);
yline(isoE_val, 'color', 'g', 'LineWidth', 1.);
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
if isfield(dataStr, 'kx');      xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
else;                           xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex');
end
%% -- Plotting the IsoK ARPES Image
nexttile;
if isfield(dataStr, 'kx')
    XM = squeeze(ZZ(:,indxX,:));
    YM = squeeze(YY(:,indxX,:));
    Data = squeeze(D1(:,indxX,:));
    [XR,YR,DataR]=Remap(XM,YM,Data);
    surf(XR,YR,DataR); view(2); shading flat; clim([min(D1(:)) max(D1(:))]);
elseif isfield(dataStr, 'eb') || isfield(dataStr, 'data')
    XM = squeeze(ZZ);
    YM = squeeze(YY(:,indxX,:));
    Data = squeeze(D1(:,indxX,:));
    [XR,YR,DataR]=Remap(XM,YM,Data);
    surf(XR,YR,DataR); view(2); shading flat; clim([min(D1(:)) max(D1(:))]);
else
    XR = ZZ;
    YR = YY;
    DataR = squeeze(D1(:,indxX,:));
    imagesc('XData', XR, 'YData', YR, 'CData', DataR,[min(D1(:)) max(D1(:))]);
end
img_props(0, 'Color', [1 1 1], varargin{:}); colormap(hot); 
title(sprintf('sagittal (YZ): IsoK = %.2f', isoK_val)); 
axis([min(XR(:)), max(XR(:)), min(YR(:)), max(YR(:))]);
xline(isoScan_val, 'color', 'w', 'LineWidth', 0.5, 'LineStyle','-');
yline(isoE_val, 'color', 'g', 'LineWidth', 1.);
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
%% -- Plotting the IsoE ARPES Image
nexttile;
if isfield(dataStr, 'kx')
    XM = squeeze(ZZ(indxY,:,:));
    YM = squeeze(XX(indxY,:,:));
    Data = squeeze(D1(indxY,:,:));
    [XR,YR,DataR]=Remap(XM,YM,Data);
    surf(XR,YR,DataR); view(2); shading flat; clim([min(D1(:)) max(D1(:))]);
elseif isfield(dataStr, 'eb') || isfield(dataStr, 'data')
    XM = squeeze(ZZ);
    YM = squeeze(XX(indxY,:,:));
    Data = squeeze(D1(indxY,:,:));
    [XR,YR,DataR]=Remap(XM,YM,Data);
    surf(XR,YR,DataR); view(2); shading flat; clim([min(D1(:)) max(D1(:))]);
else
    XR = ZZ;
    YR = XX;
    DataR = squeeze(D1(indxY,:,:));
    imagesc('XData', XR, 'YData', YR, 'CData', DataR,[min(D1(:)) max(D1(:))]);
end
img_props(0, 'Color', [1 1 1], varargin{:}); colormap(hot); 
title(sprintf('axial (XZ): IsoE = %.2f', isoE_val)); 
axis([min(XR(:)), max(XR(:)), min(YR(:)), max(YR(:))]);
camroll(90); ax = gca; set(ax, 'XDir', 'reverse', 'YAxisLocation', 'right'); colorbar;
xline(isoScan_val, 'color', 'w', 'LineWidth', 0.5, 'LineStyle','-');
yline(isoK_val, 'color', 'b', 'LineWidth', 1.);
if dataStr.Type == "Eb(kx,ky)"
    if isfield(dataStr, 'kx'); ylabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); xlabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex'); 
    else; ylabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex'); xlabel('$$ \bf  \tau (^{\circ}) $$', 'Interpreter', 'latex');
    end
elseif dataStr.Type == "Eb(kx,kz)"
    if isfield(dataStr, 'kx'); ylabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); xlabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex'); 
    else; ylabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex'); xlabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
    end
elseif dataStr.Type == "Eb(k,i)"
    if isfield(dataStr, 'kx'); ylabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); xlabel('$$ \bf  Index $$', 'Interpreter', 'latex'); 
    else; ylabel('$$ \bf  Index $$', 'Interpreter', 'latex'); xlabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
    end
end