function dataStr = bgrd_subtr_data(dataStr, bsub_args, plot_result)
% dataStr = bgrd_subtr_data(dataStr, bsub_args, plot_result)
%   This is a function that will background subtract / normalise the ARPES data using several
%   different approaches over a given window.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   DataC=SetContrast(Data,minFrac,maxFrac [,gamma]);
%
%   IN:
%   -   dataStr:        data structure of the ARPES data.
%   -   bsub_args:      {1×5} cell of {bgrdType,bgrdWin,bgrdScale,bgrdSmooth,normType}.
%                               bgrdType:       "AngleInt"
%                                               "AngleInt-clip"
%                                               "AngleInt-mean"
%                                               "AngleInt-pctile"
%                                               "AngleInt-mean-moving"
%                                               "AngleInt-pctile-moving"
%                                               "none"
%                               bgrdWin:        EDC angle range to
%                               integrate over to extract ISubtr.
%                               bgrdScale:   	scale coefficient to subtract (Data - intScale * ISubtr).
%                               bgrdSmooth: 	single, constant value that determine presmoothing of data to extract background. 
%                               normType:       "max"
%                                               "max-each"
%                                               "mean-each"
%                                               "none"
%   -   plot_result:  	if 1, will plot figure of the alignment, otherwise it wont.
%
%   OUT:
%   dataStr - MATLAB data structure with new additional fields below;
%   -   .meta.bsub_args:        {1×5} cell of input arguments.
%	-   .(data):                2D or 3D array of the normalised data.

%% Default parameters
if nargin < 3; plot_result = 0; end
if nargin < 2; bsub_args = {"AngleInt-clip", [], 0.25, 50, "max-each"}; end
if isempty(bsub_args); bsub_args = {"AngleInt-clip", [], 0.25, 50, "max-each"}; end
if isempty(plot_result); plot_result = 0; end
% - Extracting the fields to be used with most recent processing
[xField, yField, ~, dField] = find_data_fields(dataStr);
% - Do not compound background subtraction, so force dField to be = 'raw_data'
% dField = 'raw_data';
%% Validity checks on the input parameters
if length(bsub_args) ~= 5
    error('error: bsub_args is not a {1x5} cell - make sure you defined the arguments correctly.');
else
    bsub_args{1} = lower(string(bsub_args{1}));
    if isempty(bsub_args{2}) || ischar(bsub_args{2}) || isstring(bsub_args{2}); bsub_args{2} = mean(dataStr.(xField)(:)) + [-1,1].*0.95.*range(dataStr.(xField)(:)); end
    if isempty(bsub_args{3}) || ischar(bsub_args{3}) || isstring(bsub_args{3}); bsub_args{3} = 0.25; end
    if isempty(bsub_args{4}) || ischar(bsub_args{4}) || isstring(bsub_args{4}); bsub_args{4} = 50; end
    bsub_args{5} = char(lower(string(bsub_args{5})));
end

%% - 1 - Initialising the normalisation parameters
% - Extracting normalisation parameters
bgrdType    = bsub_args{1};
bgrdWin   	= bsub_args{2};
bgrdScale 	= bsub_args{3};
bgrdSmooth 	= bsub_args{4};
normType  	= bsub_args{5};
% - Extracting the window of integration
Win = sort([bgrdWin(1), bgrdWin(2)]);

%% - 2 - Background subtraction of the data over all scans
% - Checking the background type is valid
if strcmpi(bgrdType,"AngleInt") || strcmpi(bgrdType,"AngleInt-clip")...
        || strcmpi(bgrdType,"AngleInt-mean") || strcmpi(bgrdType,"AngleInt-pctile")...
        || strcmpi(bgrdType,"AngleInt-mean-moving") || strcmpi(bgrdType,"AngleInt-pctile-moving")...
        || strcmpi(bgrdType,"none")
else
    error('Error: bgrdType not recognised.');
end
% - Checking the normalisation type is valid
if strcmpi(normType,"max") || strcmpi(normType,"max-each") || strcmpi(normType,"mean-each") || strcmpi(normType,"none")
else; error('Error: normType not recognised.');
end
norm_data = [];
for i = 1:size(dataStr.(dField), 3)
    %% 2.1 - Executing background subtraction
    % (A) - Angle-integrated subtraction (subtracts angle-integrated spectrum)
    if strcmpi(bgrdType,"AngleInt")
        ISubtr = IntAngle(dataStr.(dField)(:,:,i), dataStr.(xField)(:,:,i), dataStr.(yField)(:,:,i), Win);
        % - HWHM adjusted to suppress uneven sensitivity of the CCD channels
        ISubtr = Gaco2(ISubtr, 0, bgrdSmooth);
        % - Subtracting background
        norm_data(:,:,i) = dataStr.(dField)(:,:,i) - bgrdScale * ISubtr;
    % (B) - Angle-integrated subtraction (anything < 0 is equal to 0)
    elseif strcmpi(bgrdType,"AngleInt-clip")
        ISubtr = IntAngle(dataStr.(dField)(:,:,i), dataStr.(xField)(:,:,i), dataStr.(yField)(:,:,i), Win);
        % - HWHM adjusted to suppress uneven sensitivity of the CCD channels
        ISubtr = Gaco2(ISubtr, 0, bgrdSmooth);
        % - Subtracting background
        norm_data(:,:,i) = dataStr.(dField)(:,:,i) - bgrdScale * ISubtr;
        % - Clipping all negativ values to zero
        norm_data(norm_data < 0) = 0;
    % (C) - Subtracting the mean value of the angle-integrated spectrum 
    elseif strcmpi(bgrdType,"AngleInt-mean")
        ISubtr = IntAngle(dataStr.(dField)(:,:,i), dataStr.(xField)(:,:,i), dataStr.(yField)(:,:,i), Win);
        % - Finding the mean value of the EDC
        ISubtr = mean(ISubtr(:));
        % - Subtracting background
        norm_data(:,:,i) = dataStr.(dField)(:,:,i) - bgrdScale * ISubtr;
        % - Clipping all negative values to zero
        norm_data(norm_data < 0) = 0;
    % (D) - Subtracting a percentile value of the angle-integrated spectrum 
    elseif strcmpi(bgrdType,"AngleInt-pctile")
        ISubtr = IntAngle(dataStr.(dField)(:,:,i), dataStr.(xField)(:,:,i), dataStr.(yField)(:,:,i), Win);
        % - Finding the mean value of the EDC
        ISubtr = mean(ISubtr(:)) + bgrdScale*range(ISubtr(:));
        % - Subtracting background
        norm_data(:,:,i) = dataStr.(dField)(:,:,i) - ISubtr;
        % - Clipping all negative values to zero
        norm_data(norm_data < 0) = 0;
    % (E) - Subtracting the mean value of all EDCs moving from left to right
    elseif strcmpi(bgrdType,"AngleInt-mean-moving")
        % - Filing through all EDCs to subtract the mean value of each EDC
        for jj = 1:size(dataStr.(dField), 2)
            ISubtr              = mean(dataStr.(dField)(:,jj,i));
            norm_data(:,jj,i)   = dataStr.(dField)(:,jj,i) - bgrdScale * ISubtr;
        end
        % - Clipping all negative evalues to zero
        norm_data(norm_data < 0) = 0;
    % (F) - Subtracting the percentile value of all EDCs moving from left to right
    elseif strcmpi(bgrdType,"AngleInt-pctile-moving")
        % - Filing through all EDCs to subtract the mean value of each EDC
        for jj = 1:size(dataStr.(dField), 2)
            ISubtr = mean(dataStr.(dField)(:,jj,i)) + bgrdScale*range(dataStr.(dField)(:,jj,i));
            norm_data(:,jj,i)   = dataStr.(dField)(:,jj,i) - bgrdScale * ISubtr;
        end
        % - Clipping all negative values to zero
        norm_data(norm_data < 0) = 0;
    % (G) - No background subtraction
    elseif bgrdType == "none"
        ISubtr = 0;
        norm_data(:,:,i) = dataStr.(dField)(:,:,i) - bgrdScale * ISubtr;
    end
    %% 2.2 - Renormalising the intensity
    if strcmpi(normType,"max-each");      norm_data(:,:,i) = norm_data(:,:,i) / max(max(norm_data(:,:,i)));
    elseif strcmpi(normType,"mean-each"); norm_data(:,:,i) = norm_data(:,:,i) / mean(mean(norm_data(:,:,i)));
    elseif strcmpi(normType,"none");      norm_data(:,:,i) = norm_data(:,:,i) / 1;
    end
end
% - Forcing maximum value across all scans to be unity if required
if strcmpi(normType,"max") || strcmpi(normType,"max-each")
    norm_data = norm_data / max(norm_data(:));
end

%% 3 - Assigning to a data structure
dataStr.data = norm_data;
dataStr.data = dataStr.data - min(dataStr.data(:));
% - Appending new data to the data object
dataStr.meta.bsub_args = bsub_args;
% - Setting NaN values to zero
dataStr.data(isnan(dataStr.data)) = 0;
% - Setting the contrast
for i = 1:size(dataStr.(dField), 3)
    dataStr.data(:,:,i) = SetContrast(dataStr.data(:,:,i),0.05, 1.);
end

%% -- For Debugging
if plot_result == 1
    fig = figure();
    fig.Position(3) = 2*450; 
    fig.Position(4) = 350;
    subplot(121); hold on;
    ImData(dataStr.raw_tht, dataStr.raw_eb, dataStr.raw_data);
    img_props();
    axis([min(dataStr.raw_tht(:)), max(dataStr.raw_tht(:)), min(dataStr.raw_eb(:)), max(dataStr.raw_eb(:))]);
    colorbar; title('Before', 'Interpreter', 'none');
    subplot(122); hold on;
    ImData(dataStr.(xField), dataStr.(yField), dataStr.data);
    img_props();
    axis([min(dataStr.(xField)(:)), max(dataStr.(xField)(:)), min(dataStr.(yField)(:)), max(dataStr.(yField)(:))]);
    colorbar; title('After', 'Interpreter', 'none');
end

end