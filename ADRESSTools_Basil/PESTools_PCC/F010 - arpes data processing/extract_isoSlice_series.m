function [sliceStrM, fig] = extract_isoSlice_series(dataStr, isoslice_args, plot_results)
% [sliceStrM, fig] = extract_isoSlice_series(dataStr, isoslice_args, plot_results)
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
%   -   isoslice_args: 	1x6 cell of {scanIndex, isoType, isoWin, isoLims, isoN, remap}.
%   -   plot_results:  	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -   cutStr:         new MATLAB data-structure.
%   -   fig:            figure output.

%% Default parameters
if nargin < 3; plot_results = 1; end
if nargin < 2; isoslice_args = cell(1,6); plot_results = 1; end 
if isempty(isoslice_args); isoslice_args = cell(1,6); end 
if isempty(plot_results); plot_results = 1; end
pp = plot_props();
fig = [];
% - Extracting the fields to be used with most recent processing
[xField, yField, ~, ~] = find_data_fields(dataStr);

%% - 1 - Initialising input parameters
% - Extracting the input parameters
scanIndex    = isoslice_args{1}; if isempty(scanIndex);    scanIndex = 1; end           % Scan index to be cut
isoType      = isoslice_args{2}; if isempty(isoType);      isoType = "IsoE"; end        % Type of slice to make: "IsoE" (horizontal) or "IsoK" (vertical)
isoWin       = isoslice_args{3}; if isempty(isoWin);       isoWin = 0.02; end           % single, constant value of the integration window.
isoLims      = isoslice_args{4}; if isempty(isoLims);      isoLims = [-1, 1]*1e5; end	% 1x2 vector which gives the window limits of the series [min, max]
isoN         = isoslice_args{5}; if isempty(isoN);         isoN = 10; end               % single, constant value that determines the total number of slices to extract
remap    	= isoslice_args{6};  if isempty(remap);        remap       = 1; end      	% Remap the slice onto a square grid? 1 or 0
% - Verifying the min/max values of the limits
isoLims = sort(isoLims);
step_size = 0.01;
if isoType == "IsoE"
    if isoLims(1) < min(dataStr.(yField)(:)) || isoLims(1) > max(dataStr.(yField)(:)); isoLims(1) = min(dataStr.(yField)(:)) + step_size; end
    if isoLims(2) < min(dataStr.(yField)(:)) || isoLims(2) > max(dataStr.(yField)(:)); isoLims(2) = max(dataStr.(yField)(:)) - step_size; end
    XX      = data_crop3D(dataStr, [], isoLims);
    max_N   = size(XX.(xField), 1);
elseif isoType == "IsoK"
    if isoLims(1) < min(dataStr.(xField)(:)) || isoLims(1) > max(dataStr.(xField)(:)); isoLims(1) = min(dataStr.(xField)(:)) + step_size; end
    if isoLims(2) < min(dataStr.(xField)(:)) || isoLims(2) > max(dataStr.(xField)(:)); isoLims(2) = max(dataStr.(xField)(:)) - step_size; end
    XX      = data_crop3D(dataStr, isoLims, []);
    max_N   = size(XX.(xField), 2);
end
% - Verifying the max value of the number of slices
if isoN > max_N; isoN = max_N; end
% - Extracting the central location for each Cut
isoVals    = linspace(isoLims(1), isoLims(2), isoN);

%% 2 - Extracting the sequence of all line profiles
for i = 1:length(isoVals)
    slice_window          = isoVals(i) + [-1, 1] * isoWin;
    [all_slices{i}, ~]    = extract_isoSlice(dataStr, {scanIndex, isoType, slice_window, remap}, plot_results);
end

%% 4 - Appending data to MATLAB data structure
sliceStrM                       = struct();
sliceStrM.isoslice_args         = isoslice_args;
sliceStrM.xField                = string(xField);
sliceStrM.FileName              = dataStr.FileName;
sliceStrM.Type                  = dataStr.Type;
sliceStrM.IsoType               = isoType;
sliceStrM.IsoVal                = isoVals;
sliceStrM.IsoRan                = isoWin;
sliceStrM.IsoNum                = isoN;
for i = 1:length(isoVals)
    sliceStrM.IsoWin                = isoVals(i) + [-1, 1] * isoWin;
    sliceStrM.XSlice{i}         	= all_slices{i}.XSlice;
    sliceStrM.YSlice{i}         	= all_slices{i}.YSlice;
    sliceStrM.DSlice{i}         	= all_slices{i}.DSlice;
end
end