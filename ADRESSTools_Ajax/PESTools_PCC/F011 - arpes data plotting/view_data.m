function fig = view_data(dataStr, isoslice_args)
% fig = view_data(dataStr, isoslice_args)
%   Use this as the general function to view all ARPES data, both during
%   and after the processing of the data. The function identifies whether
%   the APRES data is in the form of D(X,Y) for 2D data, or D(X,Y,Z) for 3D 
%   data and plots different figures accordingly. The 2D figure takes the
%   form of an image in a single figure, whereas the 3D figure 
%   is plotted has a plot of the 2D ARPES image and a slice along a 
%   particular dimension. This function ensures that the most recent 
%   processing is also shown in the figure, so you can use this before and
%   after processing to see the changes you have made to the data.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   gca_properties(type)
%   -   [h=] ImData(X,Y,Z[,style])
%   -   view_isoScan(dataStr, isoscan_args)
%   -   view_isoSlice(dataStr, isoslice_args)
%
%   IN:
%   -   dataStr:            data structure of the ARPES data.
%   -   isoslice_args:   	1x4 cell of {scanIndex, isoType, isoWin, remap}. Not required for 2D data, but can be used for 3D data.
%
%   OUT:
%   -   fig:                MATLAB figure object with the ARPES data plotted.

%% Default parameters
if nargin < 2; isoslice_args = []; end
if isempty(isoslice_args); isoslice_args = [];  end

%% 1 - Plotting ARPES data
% - For an Eb(k) scan type
if dataStr.Type == "Eb(k)"
    [~, fig] = extract_isoScan(dataStr);
% - For a 3D Eb(kx,ky) scan type
elseif dataStr.Type == "Eb(kx,ky)" || dataStr.Type == "Eb(kx,kz)" || dataStr.Type == "Eb(k,i)"
    [~, fig] = extract_isoSlice(dataStr, isoslice_args);
end
end