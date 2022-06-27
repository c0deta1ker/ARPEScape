function [cutStrM, fig] = extract_isoCut_series(dataStr, isocut_args, plot_results)
% [cutStrM, fig] = extract_isoCut_series(dataStr, isocut_args, plot_results)
%   This function plots a EDC/MDC cut through the ARPES data in 
%   D(X,Y,cutIndx) into a single figure that includes the Eb(k) 
%   dispersion. General function to view a series of cuts through ARPES
%   data.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   [h=] ImData(X,Y,Z[,style])
%   -   [XCut,DCut] = Cut(ACorr,ECorr,Data,xMode,Win)
%
%   IN:
%   -   dataStr:      	data structure of the ARPES data.
%   -   isocut_args:  	1x6 cell of {scanIndex, isocutType, isocutWin, isocutXLims, isocutYLims, isocutN}.
%   -   plot_results:  	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -   cutStr:         new MATLAB data-structure.
%   -   fig:            figure output.

%% Default parameters
if nargin < 3; plot_results = 1; end
if nargin < 2; isocut_args = cell(1,6); plot_results = 1; end 
if isempty(isocut_args); isocut_args = cell(1,6); end 
if isempty(plot_results); plot_results = 1; end
pp = plot_props();
fig = [];

%% - 1 - Initialising input parameters
% - Extracting the input parameters
scanIndex    = isocut_args{1}; if isempty(scanIndex);  scanIndex = 1; end          % Scan index to be cut
isoType      = isocut_args{2}; if isempty(isoType);    isoType = "MDC"; end  	% Type of line profile to extract. "EDC" (vertical) or "MDC" (horizontal)
isoWin       = isocut_args{3}; if isempty(isoWin);     isoWin = 0.02; end    	% single, constant value of the integration window.
isoXLims     = isocut_args{4}; if isempty(isoXLims);   isoXLims = [-1, 1]*1e5; end	% 1x2 vector which gives the x-axis limits of the train [minX, maxX]
isoYLims     = isocut_args{5}; if isempty(isoYLims);   isoYLims = [-1, 1]*1e5; end	% 1x2 vector which gives the y-axis limits of the train [minY, maxY]
isoN         = isocut_args{6}; if isempty(isoN);       isoN = 10; end        	% single, constant value that determines the total number of line profiles to extract
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
% - Verifying the min/max values of the limits
isoXLims = sort(isoXLims);
isoYLims = sort(isoYLims);
step_size = 0.01;
if isoXLims(1) < min(XDat(:)) || isoXLims(1) > max(XDat(:)); isoXLims(1) = min(XDat(:)) + step_size; end
if isoXLims(2) < min(XDat(:)) || isoXLims(2) > max(XDat(:)); isoXLims(2) = max(XDat(:)) - step_size; end
if isoYLims(1) < min(YDat(:)) || isoYLims(1) > max(YDat(:)); isoYLims(1) = min(YDat(:)) + step_size; end
if isoYLims(2) < min(YDat(:)) || isoYLims(2) > max(YDat(:)); isoYLims(2) = max(YDat(:)) - step_size; end
% - Extracting the maximum number of line-profiles possible
[xx, yy, ~] = data_crop2D(XDat, YDat, DDat, isoXLims, isoYLims);
% - Extracting the central location for each Cut
if isoType == "MDC"
    max_N = length(xx);
    if isoN > max_N; isoN = max_N; end
    isoVals    = linspace(isoYLims(1), isoYLims(2), isoN);
elseif isoType == "EDC"
    max_N = length(yy);
    if isoN > max_N; isoN = max_N; end
    isoVals    = linspace(isoXLims(1), isoXLims(2), isoN);
end

%% 2 - Extracting the sequence of all line profiles
for i = 1:length(isoVals)
    % Extracting the cut
    cut_window          = isoVals(i) + [-1, 1] * isoWin;
    [all_cuts{i}, ~]    = extract_isoCut(dataStr, {scanIndex, isoType, cut_window}, 0);
    % Cropping the cut to be consistent with limits
    if isoType == "MDC"
        [all_cuts{i}.XCut, ~, all_cuts{i}.DCut] = data_crop2D(all_cuts{i}.XCut, [], all_cuts{i}.DCut, isoXLims, []);
    elseif isoType == "EDC"
        [~, all_cuts{i}.XCut, all_cuts{i}.DCut] = data_crop2D([], all_cuts{i}.XCut, all_cuts{i}.DCut, [], isoYLims);
    end
end

%% 3.0 - Plotting the figure if required
if plot_results == 1
    fig = figure(); 
    fig.Position(3) = 2*pp.fig5x4(1); 
    fig.Position(4) = pp.fig5x4(2);
    isoCols	= jet(length(isoVals));     % Setting the color for each line profile
    %% 3.1 - Plotting the ARPES data
    subplot(1,2,1); hold on;
    % -- The plot depends on how far in the analysis 
    ImData(XDat, YDat, DDat);
    img_props([], string(xField)); 
    colbar_props([0.06 0.75 0.015 0.175], 'left');
    axis([min(XDat(:)), max(XDat(:)),min(YDat(:)), max(YDat(:))]);
    % -- Plotting the value of each line cut
    for i = 1:length(isoVals)
        if isoType == "MDC"
            patch(...
                [isoXLims(1), isoXLims(2), isoXLims(2), isoXLims(1), isoXLims(1)],...
                all_cuts{i}.IsoVal + [-all_cuts{i}.IsoWin, -all_cuts{i}.IsoWin(1), all_cuts{i}.IsoWin, all_cuts{i}.IsoWin, -all_cuts{i}.IsoWin],...
                isoCols(i,:), 'linewidth', 1, 'facealpha', 0.35, 'edgecolor', 'none');
            plot(isoXLims, [1, 1]*all_cuts{i}.IsoVal, 'k-', 'color', isoCols(i,:));
        elseif isoType == "EDC"
            patch(...
                all_cuts{i}.IsoVal + [-all_cuts{i}.IsoWin, -all_cuts{i}.IsoWin, all_cuts{i}.IsoWin, all_cuts{i}.IsoWin, -all_cuts{i}.IsoWin],...
                [isoYLims(1), isoYLims(2), isoYLims(2), isoYLims(1), isoYLims(1)],...
                isoCols(i,:), 'linewidth', 1, 'facealpha', 0.35, 'edgecolor', 'none');
            plot([1, 1]*all_cuts{i}.IsoVal, isoYLims, 'k-', 'color', isoCols(i,:)); 
        end
    end
    % -- Plotting an outline over the area where the line cuts are taken
    patch(...
        [isoXLims(1), isoXLims(2), isoXLims(2), isoXLims(1), isoXLims(1)],...
        [isoYLims(1), isoYLims(1), isoYLims(2), isoYLims(2), isoYLims(1)],...
        [0 1 0], 'linewidth', 0.5, 'facealpha', 0, 'edgecolor', [0 1 0], 'linestyle', ':');
    % -- Adding title to the figure
    title(sprintf(string(dataStr.FileName) + "; scan %i", scanIndex), 'interpreter', 'none', 'fontsize', 9);
% (A) -- If the data to be saved is an IsoSlice / IsoCut /IsoScan form
    if isfield(dataStr, 'IsoType')
        if dataStr.IsoType== "IsoE"
            if dataStr.Type == "Eb(kx,ky)"
                if dataStr.xField == "kx"; xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  \tau (^{\circ}) $$', 'Interpreter', 'latex');
                end
            elseif dataStr.Type == "Eb(kx,kz)"
                if dataStr.xField == "kx"; xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
                end
            end
        elseif dataStr.IsoType== "IsoK"
            ylabel('$$ \bf  E_B (eV) $$', 'Interpreter', 'latex'); 
            if dataStr.Type == "Eb(kx,ky)"
                if dataStr.xField == "kx"; xlabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex');
                else; xlabel('$$ \bf  \tau (^{\circ}) $$', 'Interpreter', 'latex');
                end
            elseif dataStr.Type == "Eb(kx,kz)"
                if dataStr.xField == "kx"; xlabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
                end
            end
% (B) -- Else, for the usual ARPES data structure
        else
            img_props([], string(dataStr.xField));
        end
    end
    %% 3.2 - Plotting all of the line profiles
    subplot(1,2,2); hold on;
    for i = length(isoVals):-1:1
        x = all_cuts{i}.XCut;
        d  = all_cuts{i}.DCut + all_cuts{i}.IsoVal;
        plot(x, d, 'k.-', 'color', isoCols(i,:), 'LineWidth', 1.5);
    end
    gca_props();
    if isoType == "MDC"
        xlim(isoXLims);
    elseif isoType == "EDC"
        xlim(isoYLims);
    end
    ylabel('$$ \bf  Intensity (arb.) $$', 'Interpreter', 'latex'); 
    
end

%% 4 - Appending data to MATLAB data structure
cutStrM                 	= struct();
cutStrM.isocut_args         = isocut_args;
cutStrM.xField            	= string(xField);
cutStrM.FileName         	= dataStr.FileName;
cutStrM.Type            	= dataStr.Type;
cutStrM.IsoType          	= isoType;
cutStrM.IsoVal              = isoVals;
cutStrM.IsoWin           	= [-1, 1] * isoWin;
cutStrM.IsoNum              = isoN;
for i = 1:cutStrM.IsoNum
    cutStrM.XCut{i}      	= all_cuts{i}.XCut;
    cutStrM.DCut{i}      	= all_cuts{i}.DCut;
end

end