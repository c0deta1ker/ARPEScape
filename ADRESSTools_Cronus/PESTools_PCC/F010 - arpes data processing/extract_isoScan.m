function [scanStr, fig] = extract_isoScan(dataStr, isoscan_args, plot_results)
% [scanStr, fig] = extract_isoScan(dataStr, isoscan_args, plot_results)
%   This function plots the 2D ARPES spectrum in the form of D(X,Y), along
%   with a consistent formatting.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   gca_properties(type)
%   -   [h=] ImData(X,Y,Z[,style])
%
%   IN:
%   -   dataStr:        data structure of the ARPES data.
%   -   isoscan_args: 	1x1 cell of {scanIndex}.
%   -   plot_results:  	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -   scanStr:    	new MATLAB data-structure.
%   -   fig:            figure output.

%% Default parameters
if nargin < 3; plot_results = 1; end
if nargin < 2; isoscan_args =  cell(1,1); plot_results = 1; end
if isempty(isoscan_args); isoscan_args = cell(1,1); end
if isempty(plot_results); plot_results = 1; end
pp = plot_props();
fig = [];
% - Extracting the fields to be used with most recent processing
[xField, yField, ~, dField] = find_data_fields(dataStr);

%% - 1 - Initialising input parameters
isoScan         = isoscan_args{1}; if isempty(isoScan);  isoScan   = 1; end     % Scan index to be cut

%% 1 - Extract the 2D ARPES spectra
% -- The data depends on how far in the analysis 
if isfield(dataStr, 'kx')
    XScan = dataStr.(xField)(:,:,isoScan);
    YScan = dataStr.(yField)(:,:,isoScan);
    DScan = dataStr.(dField)(:,:,isoScan);
elseif isfield(dataStr, 'data') || isfield(dataStr, 'eb')
    XScan = dataStr.(xField)(:,:,isoScan);
    YScan = dataStr.(yField)(:,:,isoScan);
    DScan = dataStr.(dField)(:,:,isoScan);
else
    XScan = dataStr.(xField);
    YScan = dataStr.(yField);
    DScan = dataStr.(dField)(:,:,isoScan);
end

%% 2.0 - Initialising the figure
if plot_results == 1
    fig = figure(); 
    fig.Position(1) = pp.figpos(1);
    fig.Position(2) = pp.figpos(2);
    fig.Position(3) = pp.fig5x4(1); 
    fig.Position(4) = pp.fig5x4(2);
    % - For an Eb(k) scan type
    if dataStr.Type == "Eb(k)"
        % -- Plotting the Eb(k) figure object
        if isfield(dataStr, 'kx'); set(fig, 'Name', 'Iso-Scan - Eb(k): A, N, K');
        elseif isfield(dataStr, 'data'); set(fig, 'Name', 'Iso-Scan - Eb(k): A, N');
        elseif isfield(dataStr, 'eb'); set(fig, 'Name', 'Iso-Scan - Eb(k): A');
        else; set(fig, 'Name', 'Iso-Scan - Eb(k): Raw');
        end
    % - For a 3D Eb(kx,ky) scan type
    elseif dataStr.Type == "Eb(kx,ky)"
        % -- Plotting the Eb(kx,ky) figure object
        if isfield(dataStr, 'kx'); set(fig, 'Name', 'Iso-Scan - Eb(kx,ky): A, N, K');
        elseif isfield(dataStr, 'data'); set(fig, 'Name', 'Iso-Scan - Eb(kx,ky): A, N');
        elseif isfield(dataStr, 'eb'); set(fig, 'Name', 'Iso-Scan - Eb(kx,ky): A');
        else; set(fig, 'Name', 'Iso-Scan - Eb(kx,ky): Raw');
        end
    % - For a 3D Eb(kx,kz) scan type
    elseif dataStr.Type == "Eb(kx,kz)"
        % -- Plotting the Eb(kx,kz) figure object
        if isfield(dataStr, 'kx'); set(fig, 'Name', 'Iso-Scan - Eb(kx,kz): A, N, K');
        elseif isfield(dataStr, 'data'); set(fig, 'Name', 'Iso-Scan - Eb(k): A, N');
        elseif isfield(dataStr, 'eb'); set(fig, 'Name', 'Iso-Scan - Eb(kx,kz): A');
        else; set(fig, 'Name', 'Iso-Scan - Eb(kx,kz): Raw');
        end
    end
    %% 2.1 - Plotting the 2D ARPES spectra
    ImData(XScan, YScan, DScan);
    img_props([], string(xField));
    colbar_props();
    % -- Adding title to the figure
    title(sprintf(string(dataStr.FileName) + "; scan %i", isoScan), 'interpreter', 'none', 'fontsize', 12);
    axis([min(XScan(:)), max(XScan(:)), min(YScan(:)), max(YScan(:))]);
end


%% 3 - Appending data to MATLAB data structure
scanStr                 = struct();
scanStr.isoscan_args    = isoscan_args;
scanStr.xField          = string(xField);
scanStr.FileName        = dataStr.FileName;
scanStr.Type            = dataStr.Type;
scanStr.IsoType         = "Scan";
scanStr.IsoVal          = isoScan;
scanStr.XScan           = XScan;
scanStr.YScan           = YScan;
scanStr.DScan           = DScan;
end