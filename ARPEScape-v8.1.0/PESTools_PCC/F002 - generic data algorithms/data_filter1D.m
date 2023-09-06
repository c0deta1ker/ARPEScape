function ydat_fil = data_filter1D(ydat, filType, filArgs, plot_result)
% ydat_fil = data_filter1D(ydat, filType, filArgs, plot_result)
%   This function applies a filtering / smoothing operation on 1D data.
%
%   IN:
%   -   ydat:           [nY×1] column vector of the output range.
%   -   filType:    	string of the type of filtering to use. Default: "savitzky-golay" ("none","gaco1","savitzky-golay","movmean","movmedian","gaussian").
%   -   filArgs:       	cell-array of the filtering / smoothing arguments:
%                            -> gaco1:          1x1     {half-width of Gaussian}
%                            -> savitzky-golay: 2x1     {Order, Window Size}
%                            -> movmean:        1x1     {Window Size}
%                            -> movmedian:      1x1     {Window Size}
%                            -> gaussian:       1x1     {half-width of Gaussian}
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   ydat_fil:       [nY×1] column vector of the filtered data.

%% Default parameters
if nargin < 2; filType = "sg"; end
if nargin < 3; filArgs = {3,9}; end
if nargin < 4; plot_result=0; end
if isempty(filType); filType = "sg"; end
if isempty(filArgs);  filArgs = {3,9}; end
if isempty(plot_result); plot_result=0;  end
%% Validity checks on the input parameters
filType = string(filType);

%% 1 - Filtering the data
if strcmpi(filType,"none") || strcmpi(filType,"")
    ydat_fil = ydat;
elseif strcmpi(filType,"gaco1") || strcmpi(filType,"g1")
    ydat_fil = Gaco1(ydat, filArgs{:});
elseif strcmpi(filType,"savitzky-golay") || strcmpi(filType,"sav-gol") || strcmpi(filType,"sg")
    ydat_fil = sgolayfilt(ydat, filArgs{:});
elseif strcmpi(filType,"movmean")
    ydat_fil = smoothdata(ydat, 'movmean', filArgs{:});
elseif strcmpi(filType,"movmedian")
    ydat_fil = smoothdata(ydat, 'movmedian', filArgs{:});
elseif strcmpi(filType,"gaussian") || strcmpi(filType,"g")
    ydat_fil = smoothdata(ydat, 'gaussian', filArgs{:});
end
%% 2 - Setting NaN values to zero
ydat_fil(isnan(ydat_fil)) = 0;

%% -- For Debugging
if plot_result == 1
    xdat = 1:length(ydat); 
    figure(); hold on;
    plot(xdat, ydat, 'k:', 'linewidth', 1);
    plot(xdat, ydat_fil, 'k-', 'linewidth', 2);
    gca_props();
    xlabel(' X ', 'fontweight', 'bold');
    ylabel(' Y ', 'fontweight', 'bold');
    axis([min(xdat(:)), max(xdat(:)), min(ydat(:)), max(ydat(:))]);
end
end