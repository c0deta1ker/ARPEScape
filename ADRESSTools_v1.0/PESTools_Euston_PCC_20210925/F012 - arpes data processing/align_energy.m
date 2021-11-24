function dataStr = align_energy(dataStr, align_args, dataStr_ref)
% dataStr = align_energy(dataStr, align_args, dataStr_ref)
%   This is a function that will align the dataStr ARPES data 
%   to the Fermi-level, VBM or a CB state. It operates using the
%   AlignEF function conditions and also allows a reference data set to be
%   used for the alignment if required.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   [EAlign, EF, Fail] = AlignEF(Data, ECorr [,eWin] [,dEWin] [,dESmooth] [,feat]);
%
%   IN:
%   -   dataStr:            data structure of the ARPES data.
%   -   align_args:         1x6 cell of {Type (either "global", "scans", "ref"), 
%                                               Scans (scan indices to be aligned), 
%                                                   eWin, dEWin, dESmooth, feat}.
%   -   dataStr_ref:     	reference data structure used for Eb alignment.
%
%   OUT:
%   dataStr - MATLAB data structure with new additional fields below;
%   -   .meta.align_args:  1x6 cell of input arguments of alignment.
%	-   .(tht):         aligned 2D or 3D array of theta.
%	-   .(eb):          aligned 2D or 3D array of energy.
%	-   .(eb_shifts): 	cell array that contains all of the Eb shifts applied to this data.

%% Default parameters
if nargin < 3; dataStr_ref = dataStr; end
if nargin < 2; align_args = {"global", [], 0, [], [], 'edge'}; end
if isempty(align_args); align_args = {"global", [], 0, [], [], 'edge'}; end
if isempty(dataStr_ref); dataStr = dataStr_ref; end
% - Extracting the fields to be used with most recent processing
[~, yField, ~, dField] = find_data_fields(dataStr);

% disp('(1) Eb alignment...')
% wbar = waitbar(0., 'Executing eb alignment...', 'Name', 'align_energy');

%% - 1 - Initialising the alignment parameters
% - Extracting alignment alignType parameters
alignType   = align_args{1};
scan_indxs	= align_args{2};
% - Extracting alignment parameters
eWin        = align_args{3}; if ischar(eWin); eWin = []; end
dEWin       = align_args{4}; if ischar(dEWin); dEWin = []; end
dESmooth    = align_args{5}; if ischar(dESmooth); dESmooth = []; end
feat        = align_args{6};
% - Defining the scan indices over the energy alignment for-loop
if alignType == "global" || alignType == "ref" || alignType == "ref (per channel)" || alignType == "fit2ref"
    scan_indxs = 1:size(dataStr.(dField), 3);
elseif alignType == "scans" || alignType == "global shift via scan"
    minVal = 1; maxVal = size(dataStr.(dField), 3);
    % - Checking max/min are not exceeded
    scan_indxs(scan_indxs < minVal) = minVal;
    scan_indxs(scan_indxs > maxVal) = maxVal;
    % - Checking that there are not duplicate numbers
    scan_indxs = unique(scan_indxs);
end
% - Defining a values of eWin and interpolating if necessary
if isempty(eWin)
    eWin = [];
elseif length(eWin) == 1
    eWin = ones(length(scan_indxs)) .* eWin;
elseif length(eWin) == 2 || length(eWin) ~= length(scan_indxs)
    eWin = linspace(eWin(1), eWin(2), length(scan_indxs));
end
%% - 2 - Binding energy alignment over all scans
eb_shifts = [];

% (A) IF ALIGNING OVER SCAN INDICES ONLY (CAN ONLY BE USED AFTER "global" or "ref"
if alignType == "scans"
    ebindx = 0;
    for i = scan_indxs
%         waitbar(i/size(scan_indxs, 2), wbar, 'Aligning ARPES data...', 'Name', 'align_energy');
        % - Eb alignment over defined scans
        % -- For a single value of eWin to be used
        if isempty(eWin) || length(eWin) == 1
            [~, eb_shift, ~] = AlignEF(dataStr.(dField)(:,:,i), dataStr.eb(:,1,i), eWin, dEWin, dESmooth, feat);
        % -- For an eWin value at each scan between eWin(1) -> eWin(2)
        else
            ebindx = ebindx + 1;
            [~, eb_shift, ~] = AlignEF(dataStr.(dField)(:,:,i), dataStr.eb(:,1,i), eWin(ebindx), dEWin, dESmooth, feat);
        end
        % -- Applying the energy shift to the data
        dataStr.eb(:,:,i) = dataStr.eb(:,:,i) - eb_shift;
        % -- Appending all the energy shifts into a matrix
        eb_shifts(i) = eb_shift;
    end
    dataStr.eb_shifts{end+1} = [scan_indxs, eb_shifts];
end

% (B) IF ALIGNING GLOBALLY TO THE ARPES DATA DIRECTLY
if alignType == "global"
    aligned_energy  = [];
    for i = scan_indxs
%         waitbar(i/size(scan_indxs, 2), wbar, 'Aligning ARPES data...', 'Name', 'align_energy');
        % - For a single value of eWin to be used
        if isempty(eWin) || length(eWin) == 1
            % -- If the alignment was previously performed, apply again
            if string(yField) == "eb"
                [aligned_energy(:,1,i), eb_shift, ~] = AlignEF(dataStr.(dField)(:,:,i), dataStr.eb(:,1,i), eWin, dEWin,dESmooth, feat);
            % -- If no alignment has been performed, do it for the first time
            elseif string(yField) == "raw_eb"
                [aligned_energy(:,1,i), eb_shift, ~] = AlignEF(dataStr.(dField)(:,:,i), dataStr.raw_eb, eWin, dEWin,dESmooth, feat);
            end
        % - For an eWin value at each scan between eWin(1) -> eWin(2)
        else
            % -- If the alignment was previously performed, apply again
            if string(yField) == "eb"
                [aligned_energy(:,1,i), eb_shift, ~] = AlignEF(dataStr.(dField)(:,:,i), dataStr.eb(:,1,i), eWin(i), dEWin,dESmooth, feat);
            % -- If no alignment has been performed, do it for the first time
            elseif string(yField) == "raw_eb"
                [aligned_energy(:,1,i), eb_shift, ~] = AlignEF(dataStr.(dField)(:,:,i), dataStr.raw_eb, eWin(i), dEWin,dESmooth, feat);
            end
        end
        % - Appending all the energy shifts into a matrix
        eb_shifts(i) = eb_shift;
    end
    dataStr.eb_shifts           = {}; 
    dataStr.eb_shifts{end+1}    = eb_shifts;
    dataStr.eb                  = aligned_energy;
end

% (C) IF ALIGNING GLOBALLY TO A REFERENCE DATA-SET
if alignType == "ref"
    dataStr.eb_ref_file     = dataStr_ref.FileName; 
    aligned_energy  = [];
    for i = scan_indxs
%         waitbar(i/size(scan_indxs, 2), wbar, 'Aligning ARPES data...', 'Name', 'align_energy');
        if i > size(dataStr_ref.raw_data, 3)
            % - For a single value of eWin to be used
            if isempty(eWin) || length(eWin) == 1
                [~, ref_eb_shift, ~] = AlignEF(dataStr_ref.raw_data(:,:,end), dataStr_ref.raw_eb, eWin, dEWin, dESmooth, feat);
            % - For an eWin value at each scan between eWin(1) -> eWin(2)
            else
                [~, ref_eb_shift, ~] = AlignEF(dataStr_ref.raw_data(:,:,end), dataStr_ref.raw_eb, eWin(i), dEWin, dESmooth, feat);
            end
        else
            % - For a single value of eWin to be used
            if isempty(eWin) || length(eWin) == 1
                [~, ref_eb_shift, ~] = AlignEF(dataStr_ref.raw_data(:,:,i), dataStr_ref.raw_eb, eWin, dEWin, dESmooth, feat);
            % - For an eWin value at each scan between eWin(1) -> eWin(2)
            else
                [~, ref_eb_shift, ~] = AlignEF(dataStr_ref.raw_data(:,:,i), dataStr_ref.raw_eb, eWin(i), dEWin, dESmooth, feat);
            end
        end
        aligned_energy(:,1,i) = dataStr.raw_eb - ref_eb_shift;
        % - Appending all the energy shifts into a matrix
        eb_shifts(i) = ref_eb_shift;
    end
    dataStr.eb_shifts           = {}; 
    dataStr.eb_shifts{end+1}    = eb_shifts;
    dataStr.eb                  = aligned_energy;
end

% (D) IF ALIGNING GLOBALLY TO A REFERENCE DATA-SET, TO EACH INDIVIDUAL CHANNEL
if alignType == "ref (per channel)"
    dataStr.eb_ref_file     = dataStr_ref.FileName; 
    aligned_energy  = [];
    ref_eb_shift    = [];
    for i = scan_indxs
%         waitbar(i/size(scan_indxs, 2), wbar, 'Aligning ARPES data...', 'Name', 'align_energy');
        for ikx = 1:size(dataStr_ref.raw_data, 2)
            if i > size(dataStr_ref.raw_data, 3)
                % - For a single value of eWin to be used
                if isempty(eWin) || length(eWin) == 1
                    [~, ref_eb_shift(ikx), ~] = AlignEF(dataStr_ref.raw_data(:,ikx,end), dataStr_ref.raw_eb, eWin, dEWin, dESmooth, feat);
                % - For an eWin value at each scan between eWin(1) -> eWin(2)
                else
                    [~, ref_eb_shift(ikx), ~] = AlignEF(dataStr_ref.raw_data(:,ikx,end), dataStr_ref.raw_eb, eWin(i), dEWin, dESmooth, feat);
                end
            else
                % - For a single value of eWin to be used
                if isempty(eWin) || length(eWin) == 1
                    [~, ref_eb_shift(ikx), ~] = AlignEF(dataStr_ref.raw_data(:,ikx,i), dataStr_ref.raw_eb, eWin, dEWin, dESmooth, feat);
                % - For an eWin value at each scan between eWin(1) -> eWin(2)
                else
                    [~, ref_eb_shift(ikx), ~] = AlignEF(dataStr_ref.raw_data(:,ikx,i), dataStr_ref.raw_eb, eWin(i), dEWin, dESmooth, feat);
                end
            end
        end
        aligned_energy(:,:,i) = dataStr.raw_eb - ref_eb_shift;
        % - Appending all the energy shifts into a matrix
        eb_shifts(i,:) = ref_eb_shift;
    end
    dataStr.eb_shifts           = {}; 
    dataStr.eb_shifts{end+1}    = eb_shifts;
    dataStr.eb                  = aligned_energy;
end

% (D) IF ALIGNING GLOBALLY TO A REFERENCE DATA-SET, TO EACH INDIVIDUAL CHANNEL
if alignType == "fit2ref"
    dataStr.eb_ref_file     = dataStr_ref.FileName; 
    aligned_energy  = [];
    % -- NEED TO FIX WHEN SELF REFERENCING, AS IT ONLY USES RAW DATA
    angleWin = 0.25*[min(dataStr_ref.raw_tht(:)), max(dataStr_ref.raw_tht(:))];   
    for i = scan_indxs
%         waitbar(i/size(scan_indxs, 2), wbar, 'Aligning ARPES data...', 'Name', 'align_energy');
        if i > size(dataStr_ref.raw_data, 3)
            % -- For a single value of eWin to be used
            if isempty(eWin) || length(eWin) == 1
                [~,ref_eb_shift] = find_Ef_ADRESS(dataStr_ref.raw_tht, dataStr_ref.raw_eb,dataStr_ref.raw_data(:,:,end),angleWin, eWin+dEWin.*[-1,1]);
            % -- For an eWin value at each scan between eWin(1) -> eWin(2)
            else
                [~,ref_eb_shift] = find_Ef_ADRESS(dataStr_ref.raw_tht, dataStr_ref.raw_eb,dataStr_ref.raw_data(:,:,end),angleWin, eWin(i)+dEWin.*[-1,1]);
            end
        else
            % -- For a single value of eWin to be used
            if isempty(eWin) || length(eWin) == 1
                [~,ref_eb_shift] = find_Ef_ADRESS(dataStr_ref.raw_tht, dataStr_ref.raw_eb, dataStr_ref.raw_data(:,:,i),angleWin, eWin+dEWin.*[-1,1]);
                
            % -- For an eWin value at each scan between eWin(1) -> eWin(2)
            else
                [~,ref_eb_shift] = find_Ef_ADRESS(dataStr_ref.raw_tht, dataStr_ref.raw_eb, dataStr_ref.raw_data(:,:,i),angleWin, eWin(i)+dEWin.*[-1,1]);
            end
        end
        if isfield(dataStr, 'eb');  aligned_energy(:,:,i) = dataStr.(yField)(:,:,i) - ref_eb_shift;
        else;                       aligned_energy(:,1,i) = dataStr.(yField) - ref_eb_shift;
        end
        % - Appending all the energy shifts into a matrix
        eb_shifts(i) = ref_eb_shift;
    end
    
    
    
    if ~isfield(dataStr, 'eb_shifts');  dataStr.eb_shifts = {}; end
    dataStr.eb_shifts{end+1}    = eb_shifts;
    dataStr.eb                  = aligned_energy;
end

% (F) IF ALIGNING TO A SCAN INDEX AND THEN SHIFTING GLOBALLY BY A CONSTANT
if alignType == "global shift via scan"
    % - Extracting the scan to make the alignment to
    ebindx = scan_indxs(1);
%     waitbar(0.5, wbar, 'Aligning ARPES data...', 'Name', 'align_energy');
    % - Eb alignment over a single defined scan
    [~, eb_shift, ~] = AlignEF(dataStr.(dField)(:,:,ebindx), dataStr.eb(:,1,ebindx), eWin(1), dEWin, dESmooth, feat);
    % -- Applying the energy shift to the data globally
    dataStr.eb = dataStr.eb - eb_shift;
    % -- Appending all the energy shifts into a cell
    dataStr.eb_shifts{end+1} = [ebindx, eb_shift];
end

% - Appending new data to the data object
dataStr.meta.align_args = align_args;

%% - 3 - Remapping matrices into 2D or 3D consistent forms
% - 3.1 Remapping to 2D
% - Remapping the dataStr.eb domain to be 2D
if size(dataStr.eb,1) == 1 || size(dataStr.eb,2) == 1
   dataStr.eb=repmat(dataStr.eb,[1, size(dataStr.(dField),2)]);
end

% - 3.2 Remapping to 3D
% - Remapping the dataStr.eb domain to be 3D
if size(dataStr.eb, 3) == 1
   dataStr.eb = repmat(dataStr.eb,[1, 1, size(dataStr.(dField),3)]);
end

% - 3.3 Remap the tht matrix
if ~isfield(dataStr, 'tht')
    dataStr.tht = dataStr.raw_tht;
    % - 3.1 Remapping to 2D
    % -- Remapping the dataStr.tht / K domain to be 2D
    if size(dataStr.tht,1) == 1 || size(dataStr.tht,2) == 1
       dataStr.tht=repmat(dataStr.tht,[size(dataStr.(dField),1), 1]);
    end
    % - 3.2 Remapping to 3D
    % -- Remapping the dataStr.tht / K domain to be 3D
    if size(dataStr.tht, 3) == 1
       dataStr.tht = repmat(dataStr.tht,[1, 1, size(dataStr.(dField),3)]);
    end
elseif size(dataStr.tht, 1) == 1
    dataStr.tht = dataStr.raw_tht;
    % - 3.1 Remapping to 2D
    % -- Remapping the dataStr.tht / K domain to be 2D
    if size(dataStr.tht,1) == 1 || size(dataStr.tht,2) == 1
       dataStr.tht=repmat(dataStr.tht,[size(dataStr.(dField),1), 1]);
    end
    % - 3.2 Remapping to 3D
    % -- Remapping the dataStr.tht / K domain to be 3D
    if size(dataStr.tht, 3) == 1
       dataStr.tht = repmat(dataStr.tht,[1, 1, size(dataStr.(dField),3)]);
    end
end

%% Close wait-bar
% close(wbar);
end