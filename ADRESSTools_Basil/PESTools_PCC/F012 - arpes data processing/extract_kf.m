function kfStr = extract_kf(XCut, DCut, kF_args, plot_results)
% kfStr = extract_kf(XCut, DCut, kF_args, plot_results)
%   This function determines the value of KF from an MDC cut using the
%   gradient method.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   [h=] ImData(X,Y,Z[,style])
%
%   IN:
%   -   XCut:           1xN vector of the k domain.
%   -   DCut:           1xN vector of the spectral intensity.
%   -   kF_args:        1x3 cell of {kfType("none","spline","gauss","smooth","xps_solver"), dkSmoothType, dkSmoothVal}.
%   -   plot_results:  	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -   cutStr:         new MATLAB data-structure.
%   -   fig:            figure output.

%% Default parameters
if nargin < 4; plot_results = 1; end
if nargin < 3; kF_args = cell(1,3); plot_results = 1; end 
if isempty(kF_args); kF_args = cell(1,3); end 
if isempty(plot_results); plot_results = 1; end
pp = plot_props();
              
%% - 1 - Initialising input parameters
kfType          = kF_args{1}; if isempty(kfType); kfType = "gradient"; end            	% kf extraction method
dkSmoothType 	= kF_args{2}; if isempty(dkSmoothType); dkSmoothType = "none"; end    	% MDC pre-smoothing method
dkSmoothVal  	= kF_args{3}; if isempty(dkSmoothVal); dkSmoothVal = 1.0; end           % MDC pre-smoothing parameter

%% - 2 - Preparing the MDC cuts
% - 2.1 - Extracting raw MDC cut
XCut0    = XCut;
DCut0    = DCut;
% - 2.2 - Smoothing the MDC cut
XCut1    = XCut;
if lower(dkSmoothType) == "none"
    DCut1 = DCut0;
elseif lower(dkSmoothType) == "gaco2"
    DCut1 = Gaco2(DCut0, dkSmoothVal, dkSmoothVal);
elseif lower(dkSmoothType) == "binomial"
    coeff = ones(1, 24);
    DCut1 = filter(coeff, dkSmoothVal, DCut0);
elseif lower(dkSmoothType) == "savitzky-golay"
    DCut1 = sgolayfilt(DCut0, dkSmoothVal, 11);
end
% - 2.3 - Differentiating the MDC cut to find Kf inflexion points
% -- Finding first derivative
dx = abs(XCut0(1) - XCut0(2));
DCut2 = diff(DCut1(:), 1) ./ dx;
XCut2 = 0.5*dx+XCut0(1:length(DCut2));
if lower(dkSmoothType) == "None"
    DCut2 = DCut2;
elseif lower(dkSmoothType) == "gaco2"
    DCut2 = Gaco2(DCut2, dkSmoothVal, dkSmoothVal);
elseif lower(dkSmoothType) == "binomial"
    coeff = ones(1, 24);
    DCut2 = filter(coeff, dkSmoothVal, DCut2);
elseif lower(dkSmoothType) == "savitzky-golay"
    DCut2 = sgolayfilt(DCut2, dkSmoothVal, 11);
end

%% - 3 - Extracting kF from the MDC cut
% Extracting the index of the maxima
[~, ilhs] = max(DCut2(:));
[~, irhs] = min(DCut2(:));
% Extracting the first estimate of kF
kFlhs = XCut2(ilhs);
kFrhs = XCut2(irhs);
% Extracting a better estimate of kF
[kFlhs, ~] = find_peak_loc(XCut2, abs(DCut2), kFlhs+[-0.02, 0.02], kfType, plot_results);
[kFrhs, ~] = find_peak_loc(XCut2, abs(DCut2), kFrhs+[-0.02, 0.02], kfType, plot_results);

%% - 4 - Plotting the result if required
if plot_results == 1
    % -- Initialising the figure
    fig = figure(); set(fig, 'Name', 'kF Extraction');
    fig.Position(3) = pp.fig16x9(1); 
    fig.Position(4) = pp.fig16x9(2);
    hold on;
    % 4.1 - Plotting the MDC curves
    % - Plotting raw MDC cut
    plot(XCut0, DCut0, '.-', 'linewidth', 4, 'color', [0.5,0.5,0.5,0.5]);
    % - Plotting smoothed MDC cut
    plot(XCut1, DCut1, '-', 'linewidth', 2, 'color', [0,0,0]);
    % - General figure formatting
    ax = gca;
    % Font properties
    ax.FontName = 'Helvetica'; ax.FontWeight = 'normal'; ax.FontSize = 15;
    % Tick properties
    ax.TickLabelInterpreter = 'latex';
    ax.XMinorTick = 'off'; ax.YMinorTick = 'off';
    ax.TickDir = 'out';
    % Ruler properties
    ax.XAxisLocation = 'bottom';            % 'bottom' | 'top' | 'origin'
    ax.YAxisLocation = 'left';                   % 'left' | 'right' | 'origin'
    % Box Styling properties
    ax.LineWidth = 1.0;
    ax.Box = 'off'; ax.Layer = 'Top';
    % Axis labels and limits
    xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  MDC $$ $$\bf intensity$$ $$\bf (arb.)$$', 'Interpreter', 'latex');
    xlim([min(XCut0(:)), max(XCut0(:))]);
    ylim([0, 1.25*max(DCut0(:))]);
    title("kF Extraction", 'interpreter', 'none', 'fontsize', 12);
    
    % 4.2 - Plotting the kf positions
    line([1 1]*kFlhs, [0, 1.25*max(DCut0(:))], 'Color', [1 0 0], 'LineWidth', 2.5, 'Linestyle', '-');
    line([1 1]*kFrhs, [0, 1.25*max(DCut0(:))], 'Color', [1 0 0], 'LineWidth', 2.5, 'Linestyle', '-');
    
    % 4.3 - Plotting the MDC derivatives
    yyaxis right;
    y2col = [0.91 0.41 0.17 0.8];
    ax.YColor = y2col;
    % - Plotting the derivative of the two
    plot(XCut2, DCut2, '-', 'linewidth', 2, 'color', y2col);
    % - General formatting
    ylim([min(DCut2(:)), max(DCut2(:))]);
    ylabel('$$ \bf  dy/dx (arb.) $$', 'Interpreter', 'latex');
    % -- - Adding a legend to the figure
    h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'.-', 'linewidth', 5, 'color', [0.5,0.5,0.5,0.5]);
    h(2) = plot(NaN,NaN,'-', 'linewidth', 5, 'color', [0,0,0]);
    h(3) = plot(NaN,NaN,'-', 'linewidth', 5, 'color', [0.91 0.41 0.17]);
    legend(h, {'$$ MDC_{data} $$','$$ MDC_{smooth} $$','$$ MDC_{kF}$$ '}, 'interpreter', 'latex');
end

%% 5 - Appending data to MATLAB data structure
kfStr              = struct();
kfStr.kF_args      = kF_args;
% -- Raw MDC cut data
kfStr.XCut0        = XCut0;
kfStr.DCut0        = DCut0;
% -- Smoothed MDC cut data
kfStr.XCut1        = XCut1;
kfStr.DCut1        = DCut1;
% -- Final MDC cut data used to find kf
kfStr.XCut2        = XCut2;
kfStr.DCut2        = DCut2;
% -- kF properties
kfStr.kFlhs        = kFlhs;
kfStr.kFrhs        = kFrhs;
kfStr.KF        = (kFrhs - kFlhs);

end