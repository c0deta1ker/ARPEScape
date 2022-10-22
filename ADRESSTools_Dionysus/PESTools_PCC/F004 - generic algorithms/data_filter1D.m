function ydat_fil = data_filter1D(ydat, filType, filVal)
% ydat_fil = data_filter1D(ydat, filType, filVal)
%   This function applies a filtering / smoothing operation on 1D data.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   ydat:           [N×1] column vector of the output range.
%   -   filType:    	string of the type of filtering to use. Default: "spline" ("none","gaco1","binomial","savitzky-golay").
%   -   filVal:       	scalar value of the filtering / smoothing strength.
%   -   plot_result:    if 1, will plot figure summary, otherwise it wont.
%
%   OUT:
%   -   ydat_fil:       [N×1] column vector of the filtered data.

%% Default parameters
if nargin < 2; filType = "savitzky-golay"; end
if nargin < 3; filVal = 5; end
if isempty(filType); filType = "savitzky-golay"; end
if isempty(filVal);  filVal = 5; end
%% Validity checks on the input parameters
filType = string(filType);
if filVal < 0; filVal = 0; end

%% 1 - Filtering the data
if strcmpi(filType,"none")
    ydat_fil = ydat;
elseif strcmpi(filType,"gaco1")
    ydat_fil = Gaco1(ydat, filVal);
elseif strcmpi(filType,"binomial") || strcmpi(filType,"bino") || strcmpi(filType,"b")
    coeff = ones(1, 24);
    ydat_fil = filter(coeff, filVal, ydat);
elseif strcmpi(filType,"savitzky-golay") || strcmpi(filType,"sav-gol") || strcmpi(filType,"sg")
    ydat_fil = sgolayfilt(ydat, filVal, 11);
end
%% 2 - Setting NaN values to zero
ydat_fil(isnan(ydat_fil)) = 0;

end