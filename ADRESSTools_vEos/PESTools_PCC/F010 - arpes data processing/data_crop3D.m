function dataStr = data_crop3D(dataStr, xField_lims, yField_lims, zField_lims)
% dataStr = crop_data(dataStr, xField_lims, yField_lims, zField_lims)
%   This function crops the ARPES x-, y- and z-independent variables over
%   a given range. The function crops the most recently processed data and
%   ensures consistency across both the independent variables and data
%   matrix.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%
%   IN:
%   -   dataStr:        loaded MATLAB data structure.
% 	-   xField_lims:    [1 x 2] row vector of xField ('raw_tht', 'tht', 'kx') crop limits.
% 	-   yField_lims:    [1 x 2] row vector of yField ('raw_eb', 'eb') crop limits.
%	-   zField_lims:    [1 x 2] row vector of zField ('hv / tltM', 'kz / ky') crop limits.
%
%   OUT:
%   -   dataStr:        modified and cropped ARPES data structure.

%% Default parameters
maxLim = 1e4;
if nargin < 2; xField_lims=[-1 1]*maxLim; yField_lims=[-1 1]*maxLim; zField_lims=[-1 1]*maxLim; end
if nargin < 3; yField_lims=[-1 1]*maxLim; zField_lims=[-1 1]*maxLim; end
if nargin < 4; zField_lims=[-1 1]*maxLim; end
if isempty(xField_lims); xField_lims=[-1 1]*maxLim;  end
if isempty(yField_lims); yField_lims=[-1 1]*maxLim;  end
if isempty(zField_lims); zField_lims=[-1 1]*maxLim;  end
if length(xField_lims) < 2; xField_lims=[-1 1]*maxLim;  end
if length(yField_lims) < 2; yField_lims=[-1 1]*maxLim;  end
if length(zField_lims) < 2; zField_lims=[-1 1]*maxLim;  end
% - Sorting the cropping limits in ascending order
xField_lims = sort(xField_lims);
yField_lims = sort(yField_lims);
zField_lims = sort(zField_lims);
% - Extracting the fields to be used with most recent processing
[xField, yField, zField, dField] = find_data_fields(dataStr);

% disp('Data cropping...')
% wbar = waitbar(0.1, 'Extracting cropping indices...', 'Name', 'crop_data');

%% 1 - First determine the zField indices of the crop over scan parameter
if string(zField) == "hv" || string(zField) == "tltM" || string(zField) == "index"
    [~, zIndxL] = min(abs(dataStr.(zField)(1,:) - zField_lims(1)));
    [~, zIndxU] = min(abs(dataStr.(zField)(1,:) - zField_lims(2)));
    zField_indx = [zIndxL zIndxU];
else
    thtIndx = ceil(size(dataStr.(dField), 2)/2);
    ebIndx = ceil(size(dataStr.(dField), 1)/2);
    [~, zIndxL] = min(abs(dataStr.(zField)(ebIndx,thtIndx,:) - zField_lims(1)));
    [~, zIndxU] = min(abs(dataStr.(zField)(ebIndx,thtIndx,:) - zField_lims(2)));
    zField_indx = [zIndxL zIndxU];
end
%% 2 - Next determine the xField and yField indices across each scan
for i = 1:diff(zField_indx)+1
    z_index = i  - 1 + min(zField_indx);
%     waitbar(i/diff(zField_indx)+1, wbar, 'Extracting cropping indices...', 'Name', 'crop_data');
    %% 2.1 - xField cropping indices
    if string(xField) == "raw_tht"
        [~, xIndx_L] = min(abs(dataStr.(xField)(1,:) - xField_lims(1)));
        [~, xIndx_U] = min(abs(dataStr.(xField)(1,:) - xField_lims(2)));
        xField_indx{i} = [xIndx_L xIndx_U];
    else
        ebIndx = ceil(size(dataStr.(dField), 1)/2);
        [~, xIndx_L] = min(abs(dataStr.(xField)(ebIndx,:,z_index) - xField_lims(1)));
        [~, xIndx_U] = min(abs(dataStr.(xField)(ebIndx,:,z_index) - xField_lims(2)));
        xField_indx{i} = [xIndx_L xIndx_U];
    end
    %% - 2.2 - yField cropping indices
    if string(yField) == "raw_eb"
        [~, yIndx_L] = min(abs(dataStr.(yField)(:,1) - yField_lims(1)));
        [~, yIndx_U] = min(abs(dataStr.(yField)(:,1) - yField_lims(2)));
        yField_indx{i} = [yIndx_L yIndx_U];
    else
        thtIndx = ceil(size(dataStr.(dField), 2)/2);
        [~, yIndx_L] = min(abs(dataStr.(yField)(:,thtIndx,z_index) - yField_lims(1)));
        [~, yIndx_U] = min(abs(dataStr.(yField)(:,thtIndx,z_index) - yField_lims(2)));
        yField_indx{i} = [yIndx_L yIndx_U];
    end
end
%% 3 - Validity check that the cropping indices for x- and y- are consistent
% x-field consistency checks
xField_diff = diff(cat(1, xField_indx{:}),[], 2);
for i = 1:length(xField_diff)
%     waitbar(i/length(xField_diff), wbar, 'Validity checks on x-field indices...', 'Name', 'crop_data');
    if xField_diff(i) ~= min(xField_diff)
        xField_indx{i}(1) = xField_indx{i}(1) - floor((min(xField_diff)-xField_diff(i))/2);
        xField_indx{i}(2) = xField_indx{i}(2) + ceil((min(xField_diff)-xField_diff(i))/2);   
    end
end
%y-field consistency checks
yField_diff = diff(cat(1, yField_indx{:}),[], 2);
for i = 1:length(yField_diff)
%     waitbar(i/length(yField_diff), wbar, 'Validity checks on y-field indices...', 'Name', 'crop_data');
    if yField_diff(i) ~= min(yField_diff) && yField_diff(i) > min(yField_diff)
        yField_indx{i}(2) = yField_indx{i}(2) - abs(yField_diff(i) - min(yField_diff));        
    elseif yField_diff(i) ~= min(yField_diff) && yField_diff(i)< min(yfield_diff)
        yField_indx{i}(2) = yField_indx{i}(2) + abs(yField_diff(i) - min(yField_diff));        
    end
end
%% 4 - Filing through each scan parameter and applying a subsequent crop to x and y variables
for i = 1:diff(zField_indx)+1
    z_index = i -1 + min(zField_indx);
%     waitbar(i/diff(zField_indx)+1, wbar, 'Cropping ARPES data...', 'Name', 'crop_data');
    % - Applying cropping operations
    if string(xField) == "kx"
        crp_xField{i} = dataStr.(xField)(yField_indx{i}(1):yField_indx{i}(2), xField_indx{i}(1):xField_indx{i}(2), z_index);
        crp_yField{i} = dataStr.(yField)(yField_indx{i}(1):yField_indx{i}(2), xField_indx{i}(1):xField_indx{i}(2), z_index);
        crp_dField{i} = dataStr.(dField)(yField_indx{i}(1):yField_indx{i}(2), xField_indx{i}(1):xField_indx{i}(2), z_index);
        if dataStr.Type == "Eb(k)"; crp_zField = dataStr.(zField)(1, zField_indx(1):zField_indx(2));
        elseif dataStr.Type == "Eb(kx,ky)" || dataStr.Type == "Eb(kx,kz)"  || dataStr.Type == "Eb(k,i)"; crp_zField{i} = dataStr.(zField)(yField_indx{i}(1):yField_indx{i}(2), xField_indx{i}(1):xField_indx{i}(2), z_index); 
        end
        % Cropping other fields for consistency
        crp_thtField{i} = dataStr.tht(yField_indx{i}(1):yField_indx{i}(2), xField_indx{i}(1):xField_indx{i}(2), z_index);
    elseif string(yField) == "eb"
        crp_xField{i} = dataStr.(xField)(yField_indx{i}(1):yField_indx{i}(2), xField_indx{i}(1):xField_indx{i}(2), z_index);
        crp_yField{i} = dataStr.(yField)(yField_indx{i}(1):yField_indx{i}(2), xField_indx{i}(1):xField_indx{i}(2), z_index);
        if i == diff(zField_indx)+1; crp_zField = dataStr.(zField)(1, zField_indx(1):zField_indx(2)); end
        crp_dField{i} = dataStr.(dField)(yField_indx{i}(1):yField_indx{i}(2), xField_indx{i}(1):xField_indx{i}(2), z_index);
    elseif string(xField) == "raw_tht"
        if i == diff(zField_indx)+1
            crp_xField = dataStr.(xField)(1, xField_indx{1}(1):xField_indx{i}(2));
            crp_yField = dataStr.(yField)(yField_indx{1}(1):yField_indx{i}(2), 1);
            crp_zField = dataStr.(zField)(1, zField_indx(1):zField_indx(2));
        end
        crp_dField{i} = dataStr.(dField)(yField_indx{i}(1):yField_indx{i}(2), xField_indx{i}(1):xField_indx{i}(2), z_index);
    end
end
%% 5 - Converting cell arrays into 3D matrices of data
if iscell(crp_xField); crp_xField = cat(3, crp_xField{:});
else; crp_xField = crp_xField; end
if iscell(crp_yField); crp_yField = cat(3, crp_yField{:});
else; crp_yField = crp_yField; end
if iscell(crp_zField); crp_zField = cat(3, crp_zField{:});
else; crp_zField = crp_zField; end
if iscell(crp_dField); crp_dField = cat(3, crp_dField{:});
else; crp_dField = crp_dField; end
%% 6 - Assigning the cropped variables / data to new matrices
% - Assigning the first set of variables
dataStr.(xField) = []; dataStr.(yField) = []; dataStr.(zField) = []; dataStr.(dField) = [];
dataStr.(xField) = crp_xField;
dataStr.(yField) = crp_yField;
dataStr.(zField) = crp_zField;
dataStr.(dField) = crp_dField;
if exist('crp_thtField', 'var'); crp_thtField = cat(3, crp_thtField{:}); dataStr.tht = crp_thtField; end
% - Saving the crop to the meta data
dataStr.meta.crp_lims   = [];
dataStr.meta.crp_lims   = [xField_lims, yField_lims, zField_lims];
dataStr.meta.crp_xIndx  = xField_indx;
dataStr.meta.crp_yIndx  = yField_indx;
dataStr.meta.crp_zIndx  = zField_indx;
%% Close wait-bar
% close(wbar);
end
