function dataStr = xcorr_scans(dataStr, scancorr_args)
% dataStr = xcorr_scans(dataStr, scancorr_args)
%   This function cross-correlates the scans over the defined 
%   scan index and converts it from a 3D data set, into a 2D
%   Eb(k) dispersion. The sum is over the scan parameter and works for only
%   Eb(k,i), Eb(kx,ky) or Eb(kx,kz) scans only.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   [DataXC,OffsE] = SumScanXC(Energy,Data,maxLagE[,WinE]);
%
%   IN:
%   -   dataStr:        data structure of the ARPES data.
% 	-   scancorr_args:  1x2 cell of {scan_indxs, xcorr_dE}.
%
%   OUT:
%   -   dataStr:            ARPES data structure after the cross-correlation of scans.

%% Default parameters
if nargin < 2; scancorr_args = cell(1,2); end 
if isempty(scancorr_args); scancorr_args = cell(1,2); end 

%% - 1 - Initialising input parameters
scan_indxs 	= scancorr_args{1}; if isempty(scan_indxs); scan_indxs = []; end   
xcorr_dE    = scancorr_args{2}; if isempty(xcorr_dE); xcorr_dE = 0.5; end   
% - Extracting the fields to be used with most recent processing
[xField, yField, zField, dField] = find_data_fields(dataStr);
% - Sorting the cropping limits in ascending order
if isempty(scan_indxs) 
    N = size(dataStr.(dField), 3); 
    scan_indxs = 1:N; 
else
    scan_indxs = sort(unique(scancorr_args{1}));
    if length(scan_indxs) == 1; scan_indxs = [scan_indxs(1), scan_indxs(1)]; end
end
scan_indxs = ceil(scan_indxs);
% - Finding the indices to squeeze the 3D arrays into 2D
eb_indx      = ceil(size(dataStr.(dField), 1)/2); 
tht_indx     = ceil(size(dataStr.(dField), 2)/2); 

%% 2 - EbAlign->Normalise->kConvert fields
if isfield(dataStr, 'kx')
    % - Cross-correlating the ARPES data
    Y           = squeeze(dataStr.(yField)(:,tht_indx,scan_indxs(1)));
    D           = dataStr.(dField)(:,:,scan_indxs);
    [DataXC, ~] = SumScanXC(Y, D, xcorr_dE);
    DataXC(isnan(DataXC))   = 0;
    dataStr.(dField)        = []; dataStr.(dField) = DataXC;
    % - Squeezing the variable arrays
    dataStr.(xField)        = squeeze(dataStr.(xField)(:,:,scan_indxs(1)));
    if isfield(dataStr, 'tht'); dataStr.tht = squeeze(dataStr.tht(:,:,scan_indxs(1))); end
    dataStr.(yField)        = squeeze(dataStr.(yField)(:,:,scan_indxs(1)));
    if string(zField) ~= "index"; dataStr.(zField) = squeeze(dataStr.(zField)(eb_indx,tht_indx,scan_indxs(1))); end
    % - Converting the wave-vector into a single field
    dataStr.ky              = mean(dataStr.ky(:));
    dataStr.kz              = mean(dataStr.kz(:));
%% 3 - EbAlign->Normalise fields
elseif isfield(dataStr, 'data')
    % - Cross-correlating the ARPES data
    Y           = squeeze(dataStr.(yField)(:,tht_indx,scan_indxs(1)));
    D           = dataStr.(dField)(:,:,scan_indxs);
    [DataXC, ~] = SumScanXC(Y, D, xcorr_dE);
    DataXC(isnan(DataXC)) = 0;
    dataStr.(dField) = []; dataStr.(dField) = DataXC;
    % - Squeezing the variable arrays
    dataStr.(xField) = squeeze(dataStr.(xField)(:,:,scan_indxs(1)));
    dataStr.(yField) = squeeze(dataStr.(yField)(:,:,scan_indxs(1)));
    dataStr.(zField) = dataStr.(zField)(scan_indxs(1));
%% 4 - EbAlign fields
elseif isfield(dataStr, 'eb')
    % - Cross-correlating the ARPES data
    Y           = squeeze(dataStr.(yField)(:,tht_indx,scan_indxs(1)));
    D           = dataStr.(dField)(:,:,scan_indxs);
    [DataXC, ~] = SumScanXC(Y, D, xcorr_dE);
    DataXC(isnan(DataXC)) = 0;
    dataStr.(dField) = []; dataStr.(dField) = DataXC;
    % - Squeezing the variable arrays
    dataStr.(xField) = squeeze(dataStr.(xField)(eb_indx,:,scan_indxs(1)));
    dataStr.(yField) = squeeze(dataStr.(yField)(:,tht_indx,scan_indxs(1)));
    dataStr.(zField) = dataStr.(zField)(scan_indxs(1));
%% 5 - Raw, unprocessed data fields
else
    % - Squeezing the variable arrays
    dataStr.(zField)    = dataStr.(zField)(scan_indxs(1));
    % - Cross-correlating the ARPES data
    [DataXC, ~]         = SumScanXC(dataStr.(yField), dataStr.(dField)(:,:,scan_indxs), xcorr_dE);
    DataXC(isnan(DataXC)) = 0;
    dataStr.(dField) = []; dataStr.(dField) = DataXC;
end
%% 6 - Setting the identifier to Eb(k) now after XCorr
dataStr.meta.scancorr_args = scancorr_args;
dataStr.Type    = "Eb(k)";
dataStr.tltM    = mean(dataStr.tltM);
dataStr.hv      = mean(dataStr.hv);

end
