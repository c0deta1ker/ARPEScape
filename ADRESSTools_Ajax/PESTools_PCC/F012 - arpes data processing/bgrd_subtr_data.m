function dataStr = bgrd_subtr_data(dataStr, bsub_args)
% dataStr = bgrd_subtr_data(dataStr, bsub_args)
%   This is a function that will background subtract / normalise the ARPES data using several
%   different approaches over a given window.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   DataC=SetContrast(Data,minFrac,maxFrac [,gamma]);
%
%   IN:
%   -   dataStr:        data structure of the ARPES data.
%   -   bsub_args:      1x5 cell of {bgrdType (either "angle/kx-int (fixed edc)", "angle/kx-int (fixed mean)", "angle/kx-int (moving mean)", "none"), 
%                                       bgrdWin, bgrdScale, bgrdSmooth, 
%                                           normType (either "max (all scans)", "max (each scan)", "mean (each scan)", "none")}.
%
%   OUT:
%   dataStr - MATLAB data structure with new additional fields below;
%   -   .meta.bsub_args:        1x4 cell of input arguments.
%	-   .(data):                2D or 3D array of the normalised data.

%% Default parameters
if nargin < 2; bsub_args = {"angle/kx-int (fixed mean)", [-5, 5], 0.25, 50, "none"}; end
if isempty(bsub_args); bsub_args = {"angle/kx-int (fixed mean)", [-5, 5], 0.25, 50, "none"}; end
% - Extracting the fields to be used with most recent processing
[~, yField, ~, dField] = find_data_fields(dataStr);
% - Do not compound background subtraction, so force dField to be = 'raw_data'
% dField = 'raw_data';

%% - 1 - Initialising the normalisation parameters
% - Extracting normalisation parameters
bgrdType    = bsub_args{1};
bgrdWin   	= bsub_args{2};
bgrdScale 	= bsub_args{3};
bgrdSmooth 	= bsub_args{4};
normType  	= bsub_args{5};
% - Extracting the window of integration
Win = sort([bgrdWin(1), bgrdWin(2)]);

%% - 2 - Data normalisation over all scans
for i = 1:size(dataStr.(dField), 3)
%     waitbar(i/size(dataStr.(dField), 3), wbar, 'Normalising ARPES data...', 'Name', 'normalise_data');
    %% 2.1 - Executing background subtraction
    % (A) - Angle-integrated subtraction (subtracts angle-integrated spectrum)
    if bgrdType == "angle/kx-int (fixed edc)"
        % - Extracting the EDC (vertical line profile) over a given angular-integration range of the background
        cutType         = "EDC";      	% Type of cut to make. Either "MDC" (horizontal) or "EDC" (vertical)
        cutWin          = Win;          % Integration window of the cut to be made.
        isocut_args     = {i, cutType, cutWin};
        [cutStr, ~]     = extract_isoCut(dataStr, isocut_args, 0);
        ISubtr          = cutStr.DCut;
        % - HWHM adjusted to suppress uneven sensitivity of the CCD channels
        ISubtr = Gaco2(ISubtr, 0, bgrdSmooth);
        % - Subtracting background
        norm_data(:,:,i) = dataStr.(dField)(:,:,i) - bgrdScale .* ISubtr;
        norm_data(norm_data < 0) = 0;
%         % - Plotting the background intensity
%         figure(); hold on;
%         plot(ISubtr, dataStr.(yField)(:,1,i), 'linewidth', 2);
%         ylim([min(dataStr.(yField)(:,1,i)), max(dataStr.(yField)(:,1,i))]);
%         title('EDC Background Curves');

    % (B) - Angle-integrated subtraction (subtracts angle-integrated spectrum) (per channel)
    elseif bgrdType == "angle/kx-int (fixed edc) (per channel)"
        % - Extracting the EDC (vertical line profile) over a given angular-integration range of the background
        cutType         = "EDC";      	% Type of cut to make. Either "MDC" (horizontal) or "EDC" (vertical)
        cutWin          = Win;          % Integration window of the cut to be made.
        isocut_args     = {i, cutType, cutWin};
        [cutStr, ~]     = extract_isoCut(dataStr, isocut_args, 0);
        ISubtr          = cutStr.DCut;
        % - HWHM adjusted to suppress uneven sensitivity of the CCD channels
        ISubtr = Gaco2(ISubtr, 0, bgrdSmooth);
        % - Subtracting background
        for j = 1:size(dataStr.(dField), 2)
            norm_data(:,j,i) = dataStr.(dField)(:,j,i) - bgrdScale .* ISubtr;
        end
        norm_data(norm_data < 0) = 0;
%         % - Plotting the background intensity
%         figure(); hold on;
%         plot(ISubtr, dataStr.(yField)(:,1,i), 'linewidth', 2);
%         plot(dataStr.(dField)(:,1,i), dataStr.(yField)(:,1,i), 'linewidth', 2);
%         ylim([min(dataStr.(yField)(:,1,i)), max(dataStr.(yField)(:,1,i))]);
%         title('EDC Background Curves');
        
    % (C) - Fixed mean background subtraction (subtracts a fixed mean value over a non-dispersive region)
    elseif bgrdType == "angle/kx-int (fixed mean)"
        % - Extracting the EDC (vertical line profile) over a given angular-integration range of the background
        cutType         = "EDC";      	% Type of cut to make. Either "MDC" (horizontal) or "EDC" (vertical)
        cutWin          = Win;          % Integration window of the cut to be made.
        isocut_args     = {i, cutType, cutWin};
        [cutStr, ~]     = extract_isoCut(dataStr, isocut_args, 0);
        
        ISubtr          = cutStr.DCut;
        % - Finding the mean value of the EDC
        ISubtr        	= mean(ISubtr(:));
        % - Subtracting background as the mean value of the EDC
        norm_data(:,:,i) = dataStr.(dField)(:,:,i) - bgrdScale * ISubtr;
        norm_data(norm_data < 0) = 0;
        
    % (D) - Moving mean background subtraction (subtracts a the mean value of all EDCs over whole image)
    elseif bgrdType == "angle/kx-int (moving mean)"
        % - Filing through all EDCs to subtract the mean value of each EDC
        for jj = 1:size(dataStr.(dField), 2)
            ISubtr              = mean(dataStr.(dField)(:,jj,i));
            norm_data(:,jj,i)   = dataStr.(dField)(:,jj,i) - bgrdScale * ISubtr;
        end
        norm_data(norm_data < 0) = 0;
        
    % (E) - No background subtraction
    elseif bgrdType == "none"
        ISubtr = 0;
        norm_data(:,:,i) = dataStr.(dField)(:,:,i) - bgrdScale * ISubtr;
    end
    %% 2.2 - Renormalising the intensity
    if normType == "max (each scan)"
        norm_data(:,:,i) = norm_data(:,:,i) / max(max(norm_data(:,:,i)));
    elseif normType == "mean (each scan)"
        norm_data(:,:,i) = norm_data(:,:,i) / mean(mean(norm_data(:,:,i)));
    elseif normType == "none"
        norm_data(:,:,i) = norm_data(:,:,i) / 1;
    end
end
% - Forcing maximum value across all scans to be unity if required
if normType == "max (all scans)"
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
%% Close wait-bar
% close(wbar);
end