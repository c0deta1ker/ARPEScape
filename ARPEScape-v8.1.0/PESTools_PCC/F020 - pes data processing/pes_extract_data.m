function [pesStr, fig] = pes_extract_data(arpesStr, pes_args, plot_results)
% [pesStr, fig] = pes_extract_data(arpesStr, pes_args, plot_results)
%   This function extracts an EDC cut through XPS data, allowing the
%   angle-integrated XPS spectrum to be acquired. This is then used for all
%   subsequent XPS analysis. The user can choose to make N cuts within a
%   given window if desired. This is the function used to extract ADRESS 
%   data as XPS structure within the PESTools package.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   [h=] ImData(X,Y,Z[,style])
%   -   [XCut,DCut] = Cut(ACorr,ECorr,Data,xMode,Win)
%
%   IN:
%   -   arpesStr:      	data structure of the ARPES data.
%   -   pes_args:       1x2 cell of {cutWin, cutN}; the integration window of the cut and total number of cuts.
%   -   plot_results:  	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -   pesStr:         new MATLAB data-structure for PES data, can be used with all 'pes analysis' functions.
%   -   fig:            figure output.

%% Default parameters
if nargin < 3; plot_results = 1; end
if nargin < 2; pes_args = cell(1,2); end
if isempty(pes_args) || length(pes_args) ~= 2; pes_args = cell(1,2); end 
if isempty(plot_results); plot_results = 1; end
% - Extracting the fields to be used with most recent processing
[xField, yField, ~, dField] = find_data_fields(arpesStr);

%% - 1 - Initialising input parameters
% - Extracting the input parameters
angleWin	= pes_args{1}; if isempty(angleWin); angleWin = 0.95*[min(arpesStr.(xField)(:)), max(arpesStr.(xField)(:))]; end     % Integration window of the cut to be made.
xpsN     	= pes_args{2}; if isempty(xpsN);     xpsN = 1; end        	% single, constant value that determines the total number of line profiles to extract
% - Verifying the min/max values of the limits
Win         = sort(angleWin);
step_size   = 0.05;
if Win(1) < min(arpesStr.(xField)(:)); Win(1) = min(arpesStr.(xField)(:)) + step_size; end
if Win(2) > max(arpesStr.(xField)(:)); Win(2) = max(arpesStr.(xField)(:)) - step_size; end
% - Extracting the maximum number of line-profiles possible
[~, yy, ~] = data_crop2D(arpesStr.(xField), arpesStr.(yField), arpesStr.(dField), Win, []);
max_N = length(yy); if xpsN > max_N; xpsN = max_N; end
% - Extracting the central location for each Cut to be made
cutVals   = linspace(Win(1), Win(2), xpsN+2);
if xpsN == 1; cutRange  = 0.5*range(Win);
else; cutRange  = 0.5*range(Win)./(xpsN+1); 
end

%% 2 - Extracting all of the angle-integrated cuts
for i = 1:xpsN
    cutWin(i,:)         = cutVals(i+1) + [-1, 1] * cutRange;
    [XCut{i}, DCut{i}]  = Cut(arpesStr.(xField), arpesStr.(yField), arpesStr.(dField), 'edc', cutWin(i,:));
end

%% 3 - Appending data to MATLAB data structure
if xpsN == 1
    pesStr             = arpesStr;
    pesStr.IsoType    	= "XPS";
    pesStr.pes_args   	= {angleWin, xpsN};
    pesStr.cutWin   	= cutWin(1,:);
    pesStr.cutN     	= 1;
    pesStr.xdat         = XCut{1};
    pesStr.ydat      	= DCut{1};
else
    for i = 1:xpsN
        pesStr{i}              = arpesStr;
        pesStr{i}.IsoType    	= "XPS";
        pesStr{i}.pes_args   	= {angleWin, xpsN};
        pesStr{i}.cutWin   	= cutWin(i,:);
        pesStr{i}.cutN         = i;
        pesStr{i}.xdat         = XCut{i};
        pesStr{i}.ydat      	= DCut{i};
    end
end

%% 4.0 - Plotting the figure if required
fig = [];
if plot_results == 1
    % - Initialising the plot properties
    pp  = plot_props();
    fig = figure(); 
    fig.Position(3) = 1.5*pp.fig5x4(1); 
    fig.Position(4) = pp.fig5x4(2);
    subplot(1,3,1:2); hold on;
    % - For an Eb(k) scan type
    if arpesStr.Type == "Eb(k)"
        % -- Plotting the Eb(k) figure object
        if isfield(arpesStr, 'kx'); set(fig, 'Name', 'Iso-Cut - Eb(k): A, N, K');
        elseif isfield(arpesStr, 'data'); set(fig, 'Name', 'Iso-Cut - Eb(k): A, N');
        elseif isfield(arpesStr, 'eb'); set(fig, 'Name', 'Iso-Cut - Eb(k): A');
        else; set(fig, 'Name', 'Iso-Cut - Eb(k): Raw');
        end
    % - For a 3D Eb(kx,ky) scan type
    elseif arpesStr.Type == "Eb(kx,ky)"
        % -- Plotting the Eb(kx,ky) figure object
        if isfield(arpesStr, 'kx'); set(fig, 'Name', 'Iso-Cut - Eb(kx,ky): A, N, K');
        elseif isfield(arpesStr, 'data'); set(fig, 'Name', 'Iso-Cut - Eb(kx,ky): A, N');
        elseif isfield(arpesStr, 'eb'); set(fig, 'Name', 'Iso-Cut - Eb(kx,ky): A');
        else; set(fig, 'Name', 'Iso-Cut - Eb(kx,ky): Raw');
        end
    % - For a 3D Eb(kx,kz) scan type
    elseif arpesStr.Type == "Eb(kx,kz)"
        % -- Plotting the Eb(kx,kz) figure object
        if isfield(arpesStr, 'kx'); set(fig, 'Name', 'Iso-Cut - Eb(kx,kz): A, N, K');
        elseif isfield(arpesStr, 'data'); set(fig, 'Name', 'Iso-Cut - Eb(k): A, N');
        elseif isfield(arpesStr, 'eb'); set(fig, 'Name', 'Iso-Cut - Eb(kx,kz): A');
        else; set(fig, 'Name', 'Iso-Cut - Eb(kx,kz): Raw');
        end
    end
    % - Defining the XPS plot colors
    isoCols	= jet(xpsN);
    % - Defining the axes limits
    xLims = [min(arpesStr.(xField)(:)), max(arpesStr.(xField)(:))];
    yLims = [min(arpesStr.(yField)(:)), max(arpesStr.(yField)(:))];
    %% 4.1 - Plotting the Eb(k) at the scan value
    % -- The plot depends on how far in the analysis 
    ImData(arpesStr.(xField), arpesStr.(yField), arpesStr.(dField));
    img_props(); cbar_props();
    axis([xLims, yLims]);
    % -- Plotting an outline over the integrated sliced region
    patch([Win(1), Win(1), Win(2), Win(2), Win(1)],[-1e3, 1e3, 1e3, -1e3, -1e3],...
        [0 1 0], 'linewidth', 1, 'facealpha', 0, 'edgecolor', [0 0 1]);
    % -- Plotting the value of each line cut
    for i = 1:xpsN
        patch(...
            [cutWin(i,1), cutWin(i,1), cutWin(i,2), cutWin(i,2), cutWin(i,1)],...
            [yLims(1), yLims(2), yLims(2), yLims(1), yLims(1)],...
            isoCols(i,:), 'linewidth', 1, 'facealpha', 0.35, 'edgecolor', 'none');
        plot([1, 1]*mean(cutWin(i,:)), yLims, 'k-', 'color', isoCols(i,:)); 
    end
    % -- Adding title to the figure
    title(sprintf(string(arpesStr.FileName)), 'interpreter', 'none', 'fontsize', 9);
    %% 4.2 - Plotting the Iso-Cut
    subplot(1,3,3); hold on;
    for i = 1:xpsN
        plot(DCut{i}, XCut{i}, 'k.-', 'color', [0 0 1], 'LineWidth', 1.5, 'color', isoCols(i,:));
    end
    gca_props();
    ax = gca;
    ax.XAxisLocation = 'bottom';            % 'bottom' | 'top' | 'origin'
    ax.YAxisLocation = 'right';             % 'left' | 'right' | 'origin'
    % - Axis labels and limits
    ylabel('$$ \bf  E_{B} (eV) $$', 'Interpreter', 'latex');
    xlabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
    % - Plotting the x- and y-axes
    xl = xlim; yl = [min(arpesStr.(yField)(:)), max(arpesStr.(yField)(:))];
    axis([xl(1), 1.05*xl(2), yl(1), yl(2)]);
end

end