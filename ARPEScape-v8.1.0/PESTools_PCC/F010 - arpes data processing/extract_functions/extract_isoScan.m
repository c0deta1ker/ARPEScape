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
if nargin < 2; isoscan_args =  cell(1,1); plot_results = 1; end
if nargin < 3; plot_results = 1; end
if isempty(isoscan_args); isoscan_args = cell(1,1); end
if isempty(plot_results); plot_results = 1; end
% - Extracting the fields to be used with most recent processing
[xField, yField, ~, dField] = find_data_fields(dataStr);

%% 1 - Initialising input parameters
isoScan         = isoscan_args{1}; if isempty(isoScan);  isoScan   = 1; end     % Scan index to be cut
%% 2 - Extract the 2D ARPES spectra
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
%% 4 - Plot the results
fig = [];
if plot_results == 1
    pp = plot_props();
    fig = figure(); 
    fig.Position(3) = pp.fig5x4(1); 
    fig.Position(4) = pp.fig5x4(2);
    if isfield(dataStr, 'kx');          ana_stage = "Algnd, Nrmlzed, Kconv";
    elseif isfield(dataStr, 'data');    ana_stage = "Algnd, Nrmlzed";
    elseif isfield(dataStr, 'eb');      ana_stage = "Algnd";
    else;                               ana_stage = "Raw";
    end
    set(fig, 'Name', sprintf("%s - %s: %s", scanStr.Type, scanStr.IsoType, ana_stage));
    ImData(scanStr.XScan, scanStr.YScan, scanStr.DScan);
    img_props(); cbar_props(); 
    axis([min(scanStr.XScan(:)), max(scanStr.XScan(:)), min(scanStr.YScan(:)), max(scanStr.YScan(:))]);
    title(sprintf(string(scanStr.FileName) + "; scan %i", scanStr.IsoVal), 'interpreter', 'none');
    ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
    if isfield(dataStr, 'kx');      xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
    else;                           xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex');
    end
end
end