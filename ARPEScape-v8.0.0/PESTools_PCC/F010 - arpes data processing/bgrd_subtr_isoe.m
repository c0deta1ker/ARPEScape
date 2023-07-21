function isoeStr = bgrd_subtr_isoe(isoeStr, bsub_args)
% isoeStr = bgrd_subtr_isoe(isoeStr, bsub_args)
%   This function background subtracts isoe ARPES data, given a specific
%   window [minKx, maxKx] over which to cut the isoe surface. The entire
%   isoe surface is then divided by this region, to normalised the
%   background to be uniform across all isoe scans.
%
%   REQ. FUNCTIONS:
%   -   AA = Gaco2(A,hwX,hwY [,hsX] [,hsY])
%   -   [XCut,DCut] = Cut(ACorr,ECorr,Data,xMode,Win)
%   -   DataC=SetContrast(Data,minFrac,maxFrac [,gamma]);
%
%   IN:
%   -   isoeStr:        data structure of the sliced Iso-E ARPES data.
%   -   bsub_args:      1x5 cell of {bgrdType (either "angle/kx-int (fixed edc)", "angle/kx-int (fixed mean)", "angle/kx-int (moving mean)", "none"), 
%                                       bgrdWin, bgrdScale, bgrdSmooth, 
%                                           normType (either "init max", "final max", "none")}.
%   OUT: 
%   -   isoeStr:        output data structure after transformation.

%% Default parameters
if nargin < 2; bsub_args = {"angle/kx-int (fixed mean)", [-0.65, -0.45], 0.25, 50, "none"}; end
if isempty(bsub_args); bsub_args = {"angle/kx-int (fixed mean)", [-0.65, -0.45], 0.25, 50, "none"}; end
% disp('IsoE Background Subtraction...')
% wbar = waitbar(0.5, 'Background Subtracting IsoE data...', 'Name', 'isoe_bgrd_sbtr');

%% - 1 - Initialising the normalisation parameters
% - Extracting normalisation parameters
bgrdType    = bsub_args{1};
bgrdWin   	= bsub_args{2};
bgrdScale 	= bsub_args{3};
bgrdSmooth 	= bsub_args{4};
normType  	= bsub_args{5};
% - Extracting the window of integration
Win         = sort([bgrdWin(1), bgrdWin(2)]);
% - Extracting the initial data maximum
init_max    = max(isoeStr.DSlice(:));

%% - 2 - Executing background subtraction
% - Extracting the EDC (vertical line profile) of the background
cutType         = "EDC";      	% Type of cut to make. Either "MDC" (horizontal) or "EDC" (vertical)
cutWin          = Win;          % Integration window of the cut to be made.
isocut_args     = {[], cutType, cutWin};
[cutStr, ~]     = extract_isoCut(isoeStr, isocut_args, 0);
ISubtr          = cutStr.DCut;
% (A) - Angle-integrated subtraction (subtracts angle-integrated spectrum)
if bgrdType == "angle/kx-int (fixed edc)"
    % - HWHM adjusted to suppress uneven sensitivity of the CCD channels
    ISubtr          = Gaco2(ISubtr, 0, bgrdSmooth);
    % - Subtracting background
    isoeStr.DSlice  = isoeStr.DSlice - bgrdScale .* ISubtr;
    isoeStr.DSlice(isoeStr.DSlice < 0) = 0;
    
% (B) - Fixed mean background subtraction (subtracts a fixed mean value over a non-dispersive region)
elseif bgrdType == "angle/kx-int (fixed mean)"
    % - Finding the mean value of the EDC
    ISubtr        	= mean(ISubtr(:));
    % - Subtracting background as the mean value of the EDC
    isoeStr.DSlice  = isoeStr.DSlice - bgrdScale * ISubtr;
    isoeStr.DSlice(isoeStr.DSlice < 0) = 0;

% (C) - Moving mean background subtraction (subtracts a the mean value of all EDCs over whole image)
elseif bgrdType == "angle/kx-int (moving mean)"
    % - Filing through all EDCs to subtract the mean value of each EDC
    for jj = 1:size(isoeStr.DSlice, 2)
        ISubtr                  = mean(isoeStr.DSlice(:,jj));
        isoeStr.DSlice(:,jj)    = isoeStr.DSlice(:,jj) - bgrdScale * ISubtr;
    end
    isoeStr.DSlice(isoeStr.DSlice < 0) = 0;

elseif bgrdType == "AngleInt-Norm"
    % - HWHM adjusted to suppress uneven sensitivity of the CCD channels
    ISubtr          = Gaco2(ISubtr, 0, bgrdSmooth);
    % - Subtracting background
    isoeStr.DSlice  = isoeStr.DSlice - bgrdScale .* ISubtr;
    isoeStr.DSlice(isoeStr.DSlice < 0) = 0;
    ISliceE = isoeStr.DSlice;
    ASliceE = isoeStr.XSlice;
    % - integral intensity of scans over 90% of angular range
    nA=size(ISliceE,2); n1=round(0.05*nA); n2=round(0.95*nA);
    ISliceN=ISliceE(:,n1:n2); 
    ISliceN=sum(ISliceN,2)/(n2-n1+1);
    % - normalization
    ISliceN=ISliceE./repmat(ISliceN,1,nA);
    isoeStr.DSlice = ISliceN;
    isoeStr.DSlice(isoeStr.DSlice < 0) = 0;
    
% (D) - BEST ONE, NORMALISES INCONSISTENCIES
elseif bgrdType == "Norm"
    ISliceE = isoeStr.DSlice;
    ASliceE = isoeStr.XSlice;
    % - integral intensity of scans over 90% of angular range
    nA=size(ISliceE,2); n1=round(0.05*nA); n2=round(0.95*nA);
    ISliceN=ISliceE(:,n1:n2); 
    ISliceN=sum(ISliceN,2)/(n2-n1+1);
    % - normalization
    ISliceN=ISliceE./repmat(ISliceN,1,nA);
    isoeStr.DSlice = ISliceN;
    isoeStr.DSlice(isoeStr.DSlice < 0) = 0;
% (E) - No background subtraction
elseif bgrdType == "none"
    ISubtr = 0;
    isoeStr.DSlice = isoeStr.DSlice - bgrdScale * ISubtr;
end

%% - 3 - Normalising to a given maximum
final_max   = max(isoeStr.DSlice(:));
% (A) - Normalise to the original maximum
if normType == "init max"
    isoeStr.DSlice = (isoeStr.DSlice ./ final_max) .* init_max;
% (B) - Normalise to the current maximum
elseif normType == "final max"
    isoeStr.DSlice = (isoeStr.DSlice ./ final_max);
% (C) - No background subtraction
elseif normType == "none"
    isoeStr.DSlice = isoeStr.DSlice;
end
% - Setting NaN values to zero
isoeStr.DSlice = isoeStr.DSlice - min(isoeStr.DSlice(:));
isoeStr.DSlice(isnan(isoeStr.DSlice)) = 0;
%% Close wait-bar
% close(wbar);
end
