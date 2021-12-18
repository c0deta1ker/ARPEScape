function xpsStr = normalise_xps(xpsStr, xpsnorm_args, plot_results)
% xpsStr = normalise_xps(xpsStr, xpsnorm_args, plot_results)
%   This function is used to normalise the XPS data, making the minimum
%   intensity equal zero and the maximum intensity equal one. This is not
%   necessary, but is useful when it comes to comparing qualitative changes
%   to XPS spectra. By defining aprpropriate xpsnorm_args, you can
%   normalise to a user-defined peak / position within the XPS data.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   xpsStr:      	data structure of the XPS data.
%   -   xpsnorm_args:  	1x2 cell of {eWin, dEWin}.
%   -   plot_results:  	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -   xpsStr:         data structure of the XPS data.
%   -   fig:            figure output.

%% Default parameters
if nargin < 3; plot_results = 1; end 
if nargin < 2; xpsnorm_args = cell(1,2); plot_results = 1; end 
if isempty(xpsnorm_args); xpsnorm_args = cell(1,2); end 
if isempty(plot_results); plot_results = 1; end

%% - 1 - Initialising input parameters
eWin	= xpsnorm_args{1}; if isempty(eWin); eWin = mean(xpsStr.xdat(:)); end	% approximate Eb position of a peak in the spectrum to normalise to.
dEWin 	= xpsnorm_args{2}; if isempty(dEWin); dEWin = 4; end                % Single, constant value that defines the width of the approximate Eb position

%% 2 - Normalising the XPS data
% -- Filing through all potential XPS structures in a cell array
if length(xpsStr) == 1
    % -- Shifting the data so its minimum is zero
    xpsStr.int          = xpsStr.int - min(xpsStr.int(:));
    % -- Extracting the indices over which to find the maximum
    [~, lbIndx]         = min(abs(xpsStr.xdat - (eWin - dEWin)));
    [~, ubIndx]         = min(abs(xpsStr.xdat - (eWin + dEWin)));
    Indx                = [lbIndx, ubIndx];
    % -- Finding the maximum peak value within this range
    [maxVal, maxInd]	= max(xpsStr.int(Indx(1):Indx(2)));
    peak_index = lbIndx + maxInd - 1;
    sprintf(("BE = %.2f eV, INT = %.2f"), xpsStr.xdat(peak_index), xpsStr.int(peak_index))
    % -- Assigning prior data
    xpsStr.norm.preBE      	= xpsStr.xdat;
    xpsStr.norm.preINT    	= xpsStr.int;
    xpsStr.norm.BE_peak   	= xpsStr.xdat(peak_index);
    xpsStr.norm.INT_peak  	= xpsStr.int(peak_index);
    % -- Normalising the XPS data
    xpsStr.int          = xpsStr.int ./ maxVal;
else
    for i = 1:length(xpsStr)
        % -- Shifting the data so its minimum is zero
        xpsStr{i}.int          = xpsStr{i}.int - min(xpsStr{i}.int(:));
        % -- Extracting the indices over which to find the maximum
        [~, lbIndx]         = min(abs(xpsStr{i}.xdat - (eWin - dEWin)));
        [~, ubIndx]         = min(abs(xpsStr{i}.xdat - (eWin + dEWin)));
        Indx                = [lbIndx, ubIndx];
        % -- Finding the maximum peak value within this range
        [maxVal, maxInd]	= max(xpsStr{i}.int(Indx(1):Indx(2)));
        peak_index = lbIndx + maxInd - 1;
        sprintf(("BE = %.2f eV, INT = %.2f"), xpsStr{i}.xdat(peak_index), xpsStr{i}.int(peak_index))
        % -- Assigning prior data
        xpsStr{i}.norm.preBE      	= xpsStr{i}.xdat;
        xpsStr{i}.norm.preINT    	= xpsStr{i}.int;
        xpsStr{i}.norm.BE_peak   	= xpsStr{i}.xdat(peak_index);
        xpsStr{i}.norm.INT_peak  	= xpsStr{i}.int(peak_index);
        % -- Normalising the XPS data
        xpsStr{i}.int           = xpsStr{i}.int ./ maxVal;
    end
end

%% - 3 - Plotting normalisation figure
if plot_results == 1
    % -- Initialising the figure
    pp  = plot_props();
    % -- Filing through all potential XPS structures in a cell array
    cols = jet(length(xpsStr) + 2);
    if length(xpsStr) == 1
        fig = figure(); set(fig, 'Name', 'XPS Data from ADRESS');
        fig.Position(3) = 1.25*pp.fig5x4(1); 
        fig.Position(4) = pp.fig5x4(2);
        hold on;
        % -- Plotting the XPS data after subtraction
        plot(xpsStr.xdat, xpsStr.int, 'k-', 'color', cols(1,:), 'linewidth', pp.llwidth);
        % -- Plotting a line of the normalised peak
        h = line([1, 1]*xpsStr.norm.BE_peak, [-1, 1]*1e8,...
            'color', 'r', 'linestyle', '--', 'linewidth', 1.5);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        % -- Formatting the figure
        axis([min(xpsStr.xdat(:)), max(xpsStr.xdat(:)),0, 1.25*max(xpsStr.int(:))]);
        lgnd{1} = string(xpsStr.hv) + " eV";
        title(sprintf(string(xpsStr.FileName)), 'interpreter', 'none', 'fontsize', 9);
        % - Formatting the figure
        gca_props(); grid on;
        legend(lgnd,'Location','best');
        ylabel('$$ \bf  Intensity (counts) $$', 'Interpreter', 'latex');
        xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
    else
        for i = 1:length(xpsStr)
            fig = figure(); set(fig, 'Name', 'XPS Data from ADRESS');
            fig.Position(3) = 1.25*pp.fig5x4(1); 
            fig.Position(4) = pp.fig5x4(2);
            hold on;
            % -- Plotting the XPS data after normalisation
            plot(xpsStr{i}.xdat, xpsStr{i}.int, 'k-', 'color', cols(i,:), 'linewidth', pp.llwidth);
            % -- Plotting a line of the normalised peak
            h = line([1, 1]*xpsStr{i}.norm.BE_peak, [-1, 1]*1e8,...
                'color', 'r', 'linestyle', '--', 'linewidth', 1.5);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            % -- Formatting the figure
            lgnd{1} = string(xpsStr{i}.hv) + " eV";
            axis([min(xpsStr{i}.xdat(:)), max(xpsStr{i}.xdat(:)),0, 1.25*max(xpsStr{i}.int(:))]);
            title(sprintf(string(xpsStr{i}.FileName)), 'interpreter', 'none', 'fontsize', 9);
            % - Formatting the figure
            gca_props(); grid on;
            legend(lgnd,'Location','best');
            ylabel('$$ \bf  Intensity (counts) $$', 'Interpreter', 'latex');
            xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
        end
    end
end
end