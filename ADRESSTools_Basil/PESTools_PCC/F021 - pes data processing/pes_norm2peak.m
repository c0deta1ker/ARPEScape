function ydat_norm = pes_norm2peak(xdat, ydat, norm_args, plot_results)
% ydat_norm = pes_norm2peak(xdat, ydat, norm_args, plot_results)
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
if nargin < 4; plot_results = 1; end 
if nargin < 3; norm_args = cell(1,2); plot_results = 1; end 
if isempty(norm_args); norm_args = cell(1,2); end 
if isempty(plot_results); plot_results = 1; end

%% 1 - Initialising input parameters
eWin	= norm_args{1}; if isempty(eWin); eWin = mean(xdat(:)); end                     % approximate Eb position of a peak in the spectrum to normalise to.
dEWin 	= norm_args{2}; if isempty(dEWin); dEWin = 0.5*abs(xdat(1)-xdat(end)); end  	% Single, constant value that defines the width of the approximate Eb position

%% 2 - Normalising the XPS data
% -- Shifting the data so its minimum is zero
ydat_norm           = ydat - min(ydat(:));
% -- Extracting the indices over which to find the maximum
[~, lbIndx]         = min(abs(xdat - (eWin - dEWin)));
[~, ubIndx]         = min(abs(xdat - (eWin + dEWin)));
Indx                = [lbIndx, ubIndx];
% -- Finding the maximum peak value within this range
[maxVal, maxInd]	= max(ydat_norm(Indx(1):Indx(2)));
peak_index          = lbIndx + maxInd - 1;
if plot_results == 1 
    sprintf(("xval = %.2f eV, yval = %.2f"), xdat(peak_index), ydat_norm(peak_index)) 
end
% -- Normalising the XPS data
ydat_norm           = ydat_norm ./ maxVal;


%% - 3 - Plotting normalisation figure
if plot_results == 1
    pp  = plot_props();
    fig = figure(); set(fig, 'Name', 'XPS Data from ADRESS');
    fig.Position(3) = 1.25*pp.fig5x4(1); 
    fig.Position(4) = pp.fig5x4(2);
    hold on;
    % -- Plotting the XPS data pre-subtraction
    plot(xdat, ydat, 'k-', 'linewidth', pp.llwidth);
    % -- Plotting a line of the normalised peak window
    h = line([1, 1]*eWin, [-1, 1]*5*max(ydat(:)), 'color', 'r', 'linestyle', '--', 'linewidth', 1.5);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h = line([1, 1]*eWin - dEWin, [-1, 1]*5*max(ydat(:)), 'color', 'r', 'linestyle', ':', 'linewidth', 1.5);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h = line([1, 1]*eWin + dEWin, [-1, 1]*5*max(ydat(:)), 'color', 'r', 'linestyle', ':', 'linewidth', 1.5);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % - Formatting the figure
    gca_props(); grid on;
    axis([min(xdat(:)), max(xdat(:)),0, 1.25*max(ydat(:))]);
    ylabel('$$ \bf  Initial $$', 'Interpreter', 'latex');
    xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
    % -- Plotting the XPS data after subtraction
    yyaxis right;
    plot(xdat, ydat_norm, '-', 'linewidth', pp.llwidth);
    ylabel('$$ \bf  Normalised $$', 'Interpreter', 'latex');
end

end