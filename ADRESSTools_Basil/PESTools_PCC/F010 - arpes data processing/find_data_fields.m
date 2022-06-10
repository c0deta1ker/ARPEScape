% --- General function to find the stage of the processing
function [xField, yField, zField, dField] = find_data_fields(dataStr)
% [xField, yField, zField, dField] = find_data_fields(dataStr)
%   This function determines the most recent processed 
%   fields of the ARPES data-structure so that the correct
%   post-processing can be performed in any order. The
%   fields can be used to explicitly call a field within dataStr
%   by using dataStr.(xField) for example.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   dataStr         	loaded MATLAB data structure.
%
%   OUT:
%   -   xField:            char of the most recent x-field ('raw_tht', 'tht', 'kx')
%	-   yField:            char of the most recent y-field ('raw_eb', 'eb')
%   -   zField:            char of the most recent z-field ('hv / tltM', 'kz / ky')
%	-   dField:            char of the most recent d-field ('raw_data', 'data')

% 1 - EbAlign->Normalise->kConvert fields
if isfield(dataStr, 'kx')
    xField = 'kx'; yField = 'eb'; dField = 'data';
    if dataStr.Type == "Eb(k)"; zField = 'hv';
    elseif dataStr.Type == "Eb(kx,ky)"; zField = 'ky';
    elseif dataStr.Type == "Eb(kx,kz)"; zField = 'kz';
    elseif dataStr.Type == "Eb(k,i)"; zField = 'index';
    end
% 2 - EbAlign->Normalise fields
elseif isfield(dataStr, 'data')
    xField = 'tht'; yField = 'eb'; dField = 'data';
    if dataStr.Type == "Eb(k)"; zField = 'hv';
    elseif dataStr.Type == "Eb(kx,ky)"; zField = 'tltM';
    elseif dataStr.Type == "Eb(kx,kz)"; zField = 'hv';
    elseif dataStr.Type == "Eb(k,i)"; zField = 'index';
    end
% 3 - EbAlign fields
elseif isfield(dataStr, 'eb')
    xField = 'tht'; yField = 'eb'; dField = 'raw_data';
    if dataStr.Type == "Eb(k)"; zField = 'hv';
    elseif dataStr.Type == "Eb(kx,ky)"; zField = 'tltM';
    elseif dataStr.Type == "Eb(kx,kz)"; zField = 'hv';
    elseif dataStr.Type == "Eb(k,i)"; zField = 'index';
    end
% 4 - Raw, unprocessed data fields
else
    xField = 'raw_tht'; yField = 'raw_eb'; dField = 'raw_data';
    if dataStr.Type == "Eb(k)"; zField = 'hv';
    elseif dataStr.Type == "Eb(kx,ky)"; zField = 'tltM';
    elseif dataStr.Type == "Eb(kx,kz)"; zField = 'hv';
    elseif dataStr.Type == "Eb(k,i)"; zField = 'index';
    end
end
end
