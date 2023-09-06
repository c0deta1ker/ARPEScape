function fig = view_arpes_data_2D(dataStr, varargin)
% fig = view_arpes_data_2D(dataStr, varargin)
%   Use this as the general function to view all 2D ARPES data, both during
%   and after the processing of the data. The 2D figure takes the form of 
%   a single image.This function ensures that the most recent 
%   data processing is applied.
%
%   IN:
%   -   dataStr:    data structure of the ARPES data.
%   -   varargin:   AxProps arguments: Axes properties control the appearance and behavior of an Axes object, specified as a string scalar or character vector.
%
%   OUT:
%   -   fig:        Final MATLAB figure object. 

%% Default parameters
if nargin < 2; varargin = {}; end
if isempty(varargin); varargin = {}; end
%% Default constants
pp = plot_props();
% - Extracting the fields to be used with most recent processing
[xField, yField, ~, dField] = find_data_fields(dataStr);
%% 1 - Extract the 2D ARPES spectra
% -- The data depends on how far in the analysis 
if isfield(dataStr, 'kx') || isfield(dataStr, 'data') || isfield(dataStr, 'eb')
    XScan = dataStr.(xField)(:,:,1);
    YScan = dataStr.(yField)(:,:,1);
    DScan = dataStr.(dField)(:,:,1);
else
    XScan = dataStr.(xField);
    YScan = dataStr.(yField);
    DScan = dataStr.(dField)(:,:,1);
end
%% - 2 - Plotting the final 2D Figure
fig = figure(); 
fig.Position(3) = pp.fig5x4(1); 
fig.Position(4) = pp.fig5x4(2);
% -- Set figure name
if isfield(dataStr, 'kx');          ana_stage = "Algnd, Nrmlzed, Kconv";
elseif isfield(dataStr, 'data');    ana_stage = "Algnd, Nrmlzed";
elseif isfield(dataStr, 'eb');      ana_stage = "Algnd";
else;                               ana_stage = "Raw";
end
set(fig, 'Name', sprintf("%s - %s", dataStr.Type, ana_stage));
% -- Plot the Image Data
ImData(XScan, YScan, DScan);
img_props(1, varargin{:}); cbar_props(); 
axis([min(XScan(:)), max(XScan(:)), min(YScan(:)), max(YScan(:))]);
title(sprintf(string(dataStr.FileName)), 'interpreter', 'none');
% -- Add Axes Labels
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
if isfield(dataStr, 'kx');      xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
else;                           xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex');
end