function [contStr, fig] = extract_isoConts(dataStr, isocont_args, plot_results)
% [contStr, fig] = extract_isoConts(dataStr, isocont_args, plot_results)
%   This function extracts and plots the iso-contours of a 2D ARPES image.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   [h=] ImData(X,Y,Z[,style]);
%   -   [DSlice,XSlice] = Slice(ACorr,ECorr,Data,xMode,Win) ;
%   -   [XR,YR,DataR] = Remap(XM,YM,Data);
%
%   IN:
%   -   dataStr:            data structure of the ARPES data.
%   -   isocont_args:       1x6 cell of {scanIndex, icPOI, icMinArea, icVal, icN, icSmooth}.
%   -   plot_results:    	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -   contStr:    	new MATLAB data-structure.
%   -   fig:            figure output.

%% Default parameters
if nargin < 3; plot_results = 1; end
if nargin < 2; isocont_args = cell(1,6); plot_results = 1; end
if isempty(isocont_args); isocont_args = cell(1,6); end
if isempty(plot_results); plot_results = 1; end
pp = plot_props();
fig = [];

%% - 1 - Initialising input parameters
% - Extracting the input parameters
scanIndex       = isocont_args{1}; if isempty(scanIndex);   scanIndex   = 1; end                % Scan index to be cut
icPOI           = isocont_args{2}; if isempty(icPOI);       icPOI       = [0, 15.04]; end   	% 1x2 vector of a given point of interest, aroudn which the contours are determined
icMinArea       = isocont_args{3}; if isempty(icMinArea);   icMinArea   = 0.005; end            % single, constant value that stores all contours above the minimum enclosed area
icVal           = isocont_args{4}; if isempty(icVal);       icVal       = [0.05, 0.5]; end    	% single value or 1x2 vector that gives the normalised iso-countour value, or range [minVal, maxVal]
icN             = isocont_args{5}; if isempty(icN);         icN         = 10; end               % single, constant value that determines the total number of iso-contours to extract. For single icVal, this is ignored.
icSmooth        = isocont_args{6}; if isempty(icSmooth);    icSmooth    = [1.5, 1.5]; end    	% Gaussian pre-smoothing parameters of the data prior to extracting the iso-contours
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
% - Creating a list of all contour values to file through
icVal = sort(icVal);
if length(icVal) == 1;  icRange = icVal * max(DDat(:));
else;                   icRange = linspace(icVal(1), icVal(2), icN) .* max(DDat(:));
end

%% - 2 - Extracting the Iso-Contours
DDat_raw = DDat;
% -- Gaussian pre-smooth the data if required
DDat = Gaco2(DDat, icSmooth(1), icSmooth(2));
% -- Initialising variables
figTemp = figure(); set(figTemp, 'Visible', 'off');
XCont = {};
YCont = {};
ACont = {};
% -- Filing through all potential values of the iso-contours
for i = 1:length(icRange)
    areaCont    = [];
    xCont       = {};
    yCont       = {};
    [CMatrix, ~] = contour(XDat, YDat, DDat, [1, 1]*icRange(i), 'linewidth', 2, 'linestyle', '-');
    % - Finding the areas contained within all contours
    n = 0; ii = 1;
    sz = size(CMatrix, 2);
    if ~isempty(CMatrix)
        nn(1)           = CMatrix(2,1);
        xx{1}           = CMatrix(1,2:nn(1)+1);
        yy{1}           = CMatrix(2,2:nn(1)+1);
        areaThresh(1)   = polyarea(xx{1},yy{1});
        while n+nn(ii)+ii < sz
            n               = n + nn(ii);
            ii              = ii + 1;
            nn(ii)       	= CMatrix(2,n+ii);
            xx{ii}          = CMatrix(1,n+ii+1:n+nn(ii)+ii);
            yy{ii}          = CMatrix(2,n+ii+1:n+nn(ii)+ii);
            areaThresh(ii)  = polyarea(xx{ii}, yy{ii});
        end
    else
        xx{1}       = 0;
        yy{1}       = 0;
        areaThresh  = 0;
    end
    % -- Appling the limitations from the minimum area threshold
    indx = find(areaThresh > icMinArea);
    j = 1;
    for n = indx
        xCont{j}    = xx{n};
        yCont{j}    = yy{n};
        areaCont(j) = areaThresh(n);
        j = j + 1;
    end
    areaThresh = sum(areaCont(:));

    % -- Appling the limitations from the Point Of Interest
    if ~isempty(icPOI) && length(icPOI) == 2
        POI_indx = [];
        % - Finding the index for all contours that enclose the POI
        for k = 1:length(xCont)
            xLims = [min((xCont{k}(:))), max((xCont{k}(:)))];
            yLims = [min((yCont{k}(:))), max((yCont{k}(:)))];
            if icPOI(1) < xLims(2) && icPOI(1) > xLims(1) && icPOI(2) < yLims(2) && icPOI(2) > yLims(1)
                POI_indx(end+1) = k;
            end
        end
        % - Only keep the contours that enclose the POI
        xCont       = xCont(POI_indx);
        yCont       = yCont(POI_indx);
        areaCont    = areaCont(POI_indx);
    end
    
    % -- Storing contours in a cell array
    if ~isempty(xCont) || ~isempty(yCont) || ~isempty(areaCont)
        XCont(end+1)        = xCont;
        YCont(end+1)        = yCont;
        ACont{end+1}        = areaCont;
    end
end
close(figTemp);

%% - 3 - Extracting some statistics from the iso-contours
% -- Extracting estimates for the total size of the contours
x_length = [];
y_length = [];
areas = [];
for i = 1:length(XCont)
    x_length(i) = range(XCont{i});
    y_length(i) = range(YCont{i});
    areas(i)    = ACont{i};
end
XMag        = mean(x_length);
dXMag       = 0.5 * range(x_length);
YMag        = mean(y_length);
dYMag    	= 0.5 * range(y_length);
% -- Extracting estimates for the enclosed areas
mu_area 	= mean(areas);
err_area   	= 0.5 * range(areas);
% -- Extracting estimates for the electron number density (cm^-2)
mu_ne       = 2 * (1/(2*pi))^2 * mu_area * (1e8)^2;
err_ne    	= 2 * (1/(2*pi))^2 * err_area * (1e8)^2;
% -- Defining the span of the contours as a single enclosed vector
XSpan = [XCont{1}, fliplr(XCont{end})];
YSpan = [YCont{1}, fliplr(YCont{end})];
% -- Defining the origin position of the iso-contours
X0	= mean(cell2mat(XCont));
Y0	= mean(cell2mat(YCont));

%% 4.0 - Initialising the figure
if plot_results == 1
    fig = figure(); hold on;
    fig.Position(1) = pp.figpos(1);
    fig.Position(2) = pp.figpos(2);
    fig.Position(3) = pp.fig5x4(1); 
    fig.Position(4) = pp.fig5x4(2);
    %% 4.1 - Plotting the iamge with its iso-contours
    % -- The plot depends on how far in the analysis 
    ImData(XDat, YDat, DDat_raw);
    plot(icPOI(1), icPOI(2), 'mo', 'color', 'm', 'markerfacecolor', 'm'); 
    img_props([], string(xField));
    colbar_props();
    axis([min(XDat(:)), max(XDat(:)), min(YDat(:)), max(YDat(:))]);
    % -- Adding title to the figure
    title(sprintf(string(dataStr.FileName) + "; IsoCont"), 'interpreter', 'none', 'fontsize', 9);
    axis equal;
    axis([min(XDat(:)), max(XDat(:)), min(YDat(:)), max(YDat(:))]);
    % - Initial formatting of the axis
    ax = gca;
    % -- Font properties
    ax.FontName         = 'Helvetica'; 
    ax.FontWeight       = 'normal'; 
    ax.FontSize         = 12;
    % -- Tick properties
    ax.TickLabelInterpreter = 'latex';
    ax.XMinorTick       = 'off'; 
    ax.YMinorTick       = 'off';
    ax.TickDir          = 'in';
    ax.TickLength       = [0.02 0.025];
    % -- Ruler properties
    ax.XAxisLocation    = 'bottom';            % 'bottom' | 'top' | 'origin'
    ax.YAxisLocation    = 'left';              % 'left' | 'right' | 'origin'
    % -- Box Styling properties
    ax.Color            = [0, 0, 0];
    ax.LineWidth        = 2.0;
    ax.Box              = 'on'; 
    ax.Layer            = 'Top';
    %% 4.2 - Plotting the iso-contours
    patch(XSpan, YSpan, [0 1 0], 'EdgeAlpha', 0, 'facealpha', 0.50);
    isoCols	= jet(length(XCont));     % Setting the color for each iso contour
    for i = 1:length(XCont)
        plot(XCont{i}, YCont{i}, 'color', isoCols(i,:), 'linewidth', 2, 'linestyle', '-');
        text(max(XCont{i}(:)), max(YCont{i}(:)), char(sprintf('%i', i)),...
            'horizontalalignment', 'center', 'fontsize', 14, 'color', isoCols(i,:));
    end
    % -- Plotting the total span of the iso-contours
    plot(X0+[-1,1].*0.5.*XMag(1), Y0.*[1,1], 'm-', 'linewidth', 1.5);
    plot(X0.*[1,1], Y0+[-1,1].*0.5.*YMag(1), 'm-', 'linewidth', 1.5);
    
    %% 4.3 - Formatting the axis
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
end

%% 5 - Appending data to MATLAB data structure
contStr                 = struct();
contStr.isocont_args    = isocont_args;
contStr.xField          = string(xField);
contStr.FileName        = dataStr.FileName;
contStr.Type            = dataStr.Type;
contStr.IsoType         = "IsoCont";
contStr.IsoVal          = icRange;
contStr.IsoNum          = length(XCont);
contStr.XCont           = XCont;
contStr.YCont           = YCont;
contStr.ACont           = ACont;
% Appending vectors to plot range of contours as a patch
contStr.XSpan           = XSpan;
contStr.YSpan           = YSpan;
% Appending analysis of contours
contStr.X0              = X0;
contStr.Y0              = Y0;
contStr.XMag            = [XMag, dXMag];
contStr.YMag            = [YMag, dYMag];
contStr.XYAniso      	= [XMag./YMag, abs((XMag./YMag)*sqrt((dXMag/XMag)^2 + (dYMag/YMag)^2))];
contStr.AreaMag      	= [mu_area, err_area];
contStr.Ne              = [mu_ne, err_ne];

end