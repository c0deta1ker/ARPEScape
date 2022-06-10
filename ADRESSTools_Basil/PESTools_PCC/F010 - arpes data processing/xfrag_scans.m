function dataStr_frag = xfrag_scans(dataStr, scanfrag_args)
% dataStr_frag = xfrag_scans(dataStr, scanfrag_args)
%   This function cross-fragments the scans over the defined 
%   scan index and converts it from a 3D data set, into a 2D
%   Eb(k) dispersion. This is used to treat each independent scan as a
%   seperate ARPES data set and can be used for fitting statistics.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%
%   IN:
%   -   dataStr:        data structure of the ARPES data.
% 	-   scanfrag_args:  1x1 cell of {xfrag_n}.
%
%   OUT:
%   -   dataStr_frag:   1xN cell array of all the ARPES fragmented scans

%% Default parameters
if nargin < 2; scanfrag_args = cell(1,1); end 
if isempty(scanfrag_args); scanfrag_args = cell(1,1); end 

%% 1 - Initialising input parameters
xfrag_n = scanfrag_args{1}; if isempty(xfrag_n); xfrag_n = 1; end 
% - Extracting the fields to be used with most recent processing
[~, ~, zField, dField] = find_data_fields(dataStr);

%% 2 - Fragmenting the scans
dataStr_frag = {};
if dataStr.Type == "Eb(k,i)"
    for i = 1:length(dataStr.(zField))
        dataStr_frag{i}             = dataStr;
        dataStr_frag{i}.Type     	= "Eb(k)";
        dataStr_frag{i}.(zField)	= dataStr_frag{i}.index(i);
        dataStr_frag{i}.(dField)    = squeeze(dataStr_frag{i}.(dField)(:,:,i));
    end    
elseif dataStr.Type == "Eb(kx,ky)"
    for i = 1:length(dataStr.(zField))
        dataStr_frag{i}             = dataStr;
        dataStr_frag{i}.Type     	= "Eb(k)";
        dataStr_frag{i}.(zField)  	= dataStr_frag{i}.tltM(i);
        dataStr_frag{i}.(dField)    = squeeze(dataStr_frag{i}.(dField)(:,:,i));
    end 
elseif dataStr.Type == "Eb(kx,kz)"
    for i = 1:length(dataStr.(zField))
        dataStr_frag{i}             = dataStr;
        dataStr_frag{i}.Type     	= "Eb(k)";
        dataStr_frag{i}.(zField)  	= dataStr_frag{i}.hv(i);
        dataStr_frag{i}.(dField)    = squeeze(dataStr_frag{i}.(dField)(:,:,i));
    end 
elseif dataStr.Type == "Eb(k)"
    dataStr_frag{1}                 = dataStr;
end
end
