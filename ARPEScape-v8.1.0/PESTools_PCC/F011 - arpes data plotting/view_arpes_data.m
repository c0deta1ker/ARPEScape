function fig = view_arpes_data(dataStr, slice_args, varargin)
% fig = view_arpes_data(dataStr, slice_args, varargin)
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

%% 1 - Plotting ARPES data
% - For a 2D arpes image - Eb(k) scan type
if dataStr.Type == "Eb(k)"
    fig = view_arpes_data_2D(dataStr, varargin{:});
% - For a 3D arpes image - Eb(kx,ky), Eb(kx,kz) or Eb(kx,i) scan type
elseif dataStr.Type == "Eb(kx,ky)" || dataStr.Type == "Eb(kx,kz)" || dataStr.Type == "Eb(k,i)"
    fig = view_arpes_data_3D(dataStr, slice_args, varargin{:});
end
end