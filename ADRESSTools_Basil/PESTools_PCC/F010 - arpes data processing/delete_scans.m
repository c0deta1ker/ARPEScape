function dataStr = delete_scans(dataStr, scancorr_args)
% dataStr = delete_scans(dataStr, scancorr_args)
%   This function deletes a single or linear range of scans that contain
%   anomalies. The given input is a vector of all the scan indices of the
%   ARPES scans to be deleted. Function to delete a scan - *Eb(kx,ky) or 
%   Eb(kx,kz) scans only*
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%
%   IN:
%   -   dataStr:        loaded MATLAB data structure.
% 	-   scancorr_args:  1x1 cell of {scan_indxs}.
%
%   OUT:
%   -   dataStr:                     ARPES data structure after deleting desired scans.

%% Default parameters
dataStr.meta.scancorr_args = scancorr_args;
% - Extracting the fields to be used with most recent processing
[xField, yField, zField, dField] = find_data_fields(dataStr);
% - Sorting the cropping limits in ascending order
scan_indxs = sort(scancorr_args{1});
if length(scan_indxs) == 1; scan_indxs = [scan_indxs(1), scan_indxs(1)]; end

% disp('Deleting ARPES scans...')
% wbar = waitbar(0.5, 'Deleting scans...', 'Name', 'delete_scans');

%% Validity check on input
minVal = 1;
maxVal = size(dataStr.(dField), 3);
% - Checking max/min are not exceeded
scan_indxs(scan_indxs < minVal) = minVal;
scan_indxs(scan_indxs > maxVal) = maxVal;
% - Checking that there are not duplicate numbers
scan_indxs = unique(scan_indxs);

%% - 1 - EbAlign->Normalise->kConvert fields
if isfield(dataStr, 'kx')
    dataStr.(xField)(:,:,scan_indxs)    = [];
    dataStr.tht(:,:,scan_indxs)         = [];
    dataStr.(yField)(:,:,scan_indxs)    = [];
    dataStr.(zField)(:,:,scan_indxs)    = [];
    dataStr.(dField)(:,:,scan_indxs)    = [];
%% - 2 - EbAlign->Normalise fields
elseif isfield(dataStr, 'data')
    dataStr.(xField)(:,:,scan_indxs)    = [];
    dataStr.(yField)(:,:,scan_indxs)    = [];
    dataStr.(zField)(scan_indxs)        = [];
    dataStr.(dField)(:,:,scan_indxs)    = [];
%% - 3 - EbAlign fields
elseif isfield(dataStr, 'eb')
    dataStr.(xField)(:,:,scan_indxs)    = [];
    dataStr.(yField)(:,:,scan_indxs)    = [];
    dataStr.(zField)(scan_indxs)        = [];
    dataStr.(dField)(:,:,scan_indxs)    = [];
%% - 4 - Raw, unprocessed data fields
else
    dataStr.(zField)(scan_indxs)        = [];
    dataStr.(dField)(:,:,scan_indxs)    = [];
end
%% Close wait-bar
% close(wbar);
end
