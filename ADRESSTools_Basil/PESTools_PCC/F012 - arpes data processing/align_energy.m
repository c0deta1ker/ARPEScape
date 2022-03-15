function dataStr = align_energy(dataStr, align_args, dataStr_ref, plot_result)
% dataStr = align_energy(dataStr, align_args, dataStr_ref, plot_result)
%   This is a function that will align the 'dataStr' ARPES data-structure 
%   to the Fermi-level, VBM or a CBM state. It operates using the
%   AlignEF() function conditions, but can also load in a reference data file
%   'dataStr_ref' that it will use as the energy reference.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   [EAlign, EF, Fail] = AlignEF(Data, ECorr [,eWin] [,dEWin] [,dESmooth] [,feat]);
%
%   IN:
%   -   dataStr:            data structure of the ARPES data.
%   -   align_args:         1x6 cell of {Type (either "global", "ref", "ref (per channel)", "fit2ref", "fit2ref-50%", "fit2ref-75%", "fit2ref-5%", "fit2ref-25%", "global shift via scan", "scans"), 
%                                               Scans (scan indices to be aligned), 
%                                                   eWin, dEWin, dESmooth, feat}.
%   -   dataStr_ref:     	reference data structure used for Eb alignment.
%   -   plot_results:       if 1, will plot figure of the alignment, otherwise it wont.
%
%   OUT:
%   dataStr - MATLAB data structure with new additional fields below;
%   -   .meta.align_args:  1x6 cell of input arguments of alignment.
%	-   .(tht):         aligned 2D or 3D array of theta.
%	-   .(eb):          aligned 2D or 3D array of energy.
%	-   .(eb_shifts): 	cell array that contains all of the Eb shifts applied to this data.

%% Default parameters
if nargin < 4; plot_result = 1; end
if nargin < 3; dataStr_ref = dataStr; end
if nargin < 2; align_args = {"global", [], 0, [], [], 'edge'}; end
if isempty(align_args); align_args = {"global", [], 0, [], [], 'edge'}; end
if isempty(dataStr_ref); dataStr = dataStr_ref; end
% - Extracting the fields to be used with most recent processing
[~, yField, ~, dField] = find_data_fields(dataStr);
[xField_ref, yField_ref, ~, dField_ref] = find_data_fields(dataStr_ref);

%% - 1 - Initialising the alignment parameters
% - Extracting alignment alignType parameters
alignType   = align_args{1}; alignType = lower(string(alignType));
scan_indxs	= align_args{2};
% - Extracting alignment parameters
eWin        = align_args{3}; if ischar(eWin); eWin = []; end
dEWin       = align_args{4}; if ischar(dEWin); dEWin = []; end
dESmooth    = align_args{5}; if ischar(dESmooth); dESmooth = []; end
feat        = align_args{6}; feat = char(feat);
% - Defining the scan indices over the energy alignment for-loop
if alignType ~= "scans" && alignType ~= "global shift via scan"
    scan_indxs = 1:size(dataStr.(dField), 3);
else
    minVal = 1; maxVal = size(dataStr.(dField), 3);
    % - Checking max/min are not exceeded
    scan_indxs(scan_indxs < minVal) = minVal;
    scan_indxs(scan_indxs > maxVal) = maxVal;
    % - Checking that there are not duplicate numbers
    scan_indxs = unique(scan_indxs);
end
% - Defining the values of eWin and interpolating if necessary
if isempty(eWin)
    eWin = [];
elseif length(eWin) == 1
    eWin = ones(length(scan_indxs)) .* eWin;
elseif length(eWin) == 2 || length(eWin) ~= length(scan_indxs)
    eWin = linspace(eWin(1), eWin(2), length(scan_indxs));
end

%% - 2 - Executing binding energy alignment over all scans
%% EXECUTING GLOBAL ENERGY ALIGNMENTS
if alignType == "global" || alignType == "ref" || alignType == "ref (per channel)" || alignType == "fit2ref" || alignType == "fit2ref-50%" || alignType == "fit2ref-75%" || alignType == "fit2ref-5%" || alignType == "fit2ref-25%" || alignType == "fit2ref-10%"
    eb_shifts       = []; if ~isfield(dataStr, 'eb_shifts'); dataStr.eb_shifts = {}; end
    aligned_energy	= [];
    i_eWin          = [];
    % - Iterating through all scan indices
    for i = scan_indxs
        % -- For a single value of eWin to be used
        if isempty(eWin) || length(eWin) == 1; i_eWin = eWin;
        % -- For an eWin value at each scan between eWin(1) -> eWin(2)
        else; i_eWin = eWin(i);
        end
        % -- (A) IF ALIGNING GLOBALLY TO THE ARPES DATA DIRECTLY
        if alignType == "global"
            % --- If the alignment was previously performed, apply again
            if string(yField) == "eb";          [~, eb_shift, ~] = AlignEF(dataStr.(dField)(:,:,i), dataStr.eb(:,1,i), i_eWin, dEWin, dESmooth, feat);
            % --- If no alignment has been performed, do it for the first time
            elseif string(yField) == "raw_eb";  [~, eb_shift, ~] = AlignEF(dataStr.(dField)(:,:,i), dataStr.raw_eb, i_eWin, dEWin, dESmooth, feat);
            end
        % -- (B) IF ALIGNING GLOBALLY TO A REFERENCE DATA-SET
        elseif alignType == "ref"
            if i > size(dataStr_ref.(dField_ref), 3)
                if isempty(eWin) || length(eWin) == 1;      [~, eb_shift, ~] = AlignEF(dataStr_ref.(dField_ref)(:,:,end), dataStr_ref.(yField_ref), i_eWin, dEWin, dESmooth, feat);
                else;                                       [~, eb_shift, ~] = AlignEF(dataStr_ref.(dField_ref)(:,:,end), dataStr_ref.(yField_ref), i_eWin, dEWin, dESmooth, feat);
                end
            else
                if isempty(eWin) || length(eWin) == 1;      [~, eb_shift, ~] = AlignEF(dataStr_ref.(dField_ref)(:,:,i), dataStr_ref.(yField_ref), i_eWin, dEWin, dESmooth, feat);
                else;                                       [~, eb_shift, ~] = AlignEF(dataStr_ref.(dField_ref)(:,:,i), dataStr_ref.(yField_ref), i_eWin, dEWin, dESmooth, feat);
                end
            end
        % -- (C) IF ALIGNING GLOBALLY TO A REFERENCE DATA-SET, TO EACH INDIVIDUAL X-CHANNEL
        elseif alignType == "ref (per channel)"
            % --- Iterating through all the angle (or kx) channels
            for ikx = 1:size(dataStr_ref.(dField_ref), 2)
                % ---- Executing energy alignment for each channel
                if i > size(dataStr_ref.(dField_ref), 3)
                    if isempty(eWin) || length(eWin) == 1;  [~, eb_shift(ikx), ~] = AlignEF(dataStr_ref.(dField_ref)(:,ikx,end), dataStr_ref.(yField_ref), i_eWin, dEWin, dESmooth, feat);
                    else;                                	[~, eb_shift(ikx), ~] = AlignEF(dataStr_ref.(dField_ref)(:,ikx,end), dataStr_ref.(yField_ref), i_eWin, dEWin, dESmooth, feat);
                    end
                else
                    if isempty(eWin) || length(eWin) == 1;	[~, eb_shift(ikx), ~] = AlignEF(dataStr_ref.(dField_ref)(:,ikx,i), dataStr_ref.(yField_ref), i_eWin, dEWin, dESmooth, feat);
                    else;                                 	[~, eb_shift(ikx), ~] = AlignEF(dataStr_ref.(dField_ref)(:,ikx,i), dataStr_ref.(yField_ref), i_eWin, dEWin, dESmooth, feat);
                    end
                end
            end
        % -- (D) IF ALIGNING TO THE FERMI-EDGE AND FITTING TO IT DIRECTLY
        elseif alignType == "fit2ref" || alignType == "fit2ref-50%" || alignType == "fit2ref-75%" || alignType == "fit2ref-5%" || alignType == "fit2ref-25%" || alignType == "fit2ref-10%"
            % --- Initialising variables
            if alignType == "fit2ref";          angleWin = 0.95*[min(dataStr_ref.(xField_ref)(:)), max(dataStr_ref.(xField_ref)(:))];
            elseif alignType == "fit2ref-25%";  angleWin = 0.25*[min(dataStr_ref.(xField_ref)(:)), max(dataStr_ref.(xField_ref)(:))]; 
            elseif alignType == "fit2ref-50%";  angleWin = 0.50*[min(dataStr_ref.(xField_ref)(:)), max(dataStr_ref.(xField_ref)(:))]; 
            elseif alignType == "fit2ref-75%";  angleWin = 0.75*[min(dataStr_ref.(xField_ref)(:)), max(dataStr_ref.(xField_ref)(:))];
            elseif alignType == "fit2ref-5%";   angleWin = 0.05*[min(dataStr_ref.(xField_ref)(:)), max(dataStr_ref.(xField_ref)(:))];
            elseif alignType == "fit2ref-10%";  angleWin = 0.10*[min(dataStr_ref.(xField_ref)(:)), max(dataStr_ref.(xField_ref)(:))];
            end
            % ---- Executing energy alignment by fitting to FDD
            if i > size(dataStr_ref.raw_data, 3)
                if isempty(eWin) || length(eWin) == 1;  [~,eb_shift] = find_Ef_ADRESS(dataStr_ref.(xField_ref), dataStr_ref.(yField_ref), dataStr_ref.(dField_ref)(:,:,end), angleWin, i_eWin+dEWin.*[-1,1], plot_result);
                else;                                   [~,eb_shift] = find_Ef_ADRESS(dataStr_ref.(xField_ref), dataStr_ref.(yField_ref), dataStr_ref.(dField_ref)(:,:,end), angleWin, i_eWin+dEWin.*[-1,1], plot_result);
                end
            else
                if isempty(eWin) || length(eWin) == 1;  [~,eb_shift] = find_Ef_ADRESS(dataStr_ref.(xField_ref), dataStr_ref.(yField_ref), dataStr_ref.(dField_ref)(:,:,i), angleWin, i_eWin+dEWin.*[-1,1], plot_result);
                else;                                   [~,eb_shift] = find_Ef_ADRESS(dataStr_ref.(xField_ref), dataStr_ref.(yField_ref), dataStr_ref.(dField_ref)(:,:,i), angleWin, i_eWin+dEWin.*[-1,1], plot_result);
                end
            end
        end
        %% - Appending all the aligned data vectors into a matrix
        if isfield(dataStr, 'eb');  aligned_energy(:,:,i) = dataStr.eb(:,:,i) - eb_shift;
        else;                       aligned_energy(:,1,i) = dataStr.raw_eb - eb_shift;
        end
        % - Appending all the energy shifts into a matrix
        eb_shifts(i,:) = eb_shift;
    end
    % Assigning final variables to the data structure
    dataStr.eb_ref_file         = dataStr_ref.FileName; 
    dataStr.eb_shifts{end+1}    = eb_shifts;
    dataStr.eb                  = aligned_energy;
end

%% EXECUTING SCAN SELECTED ENERGY ALIGNMENTS
if alignType == "global shift via scan" || alignType == "scans"
    if ~isfield(dataStr, 'eb'); error('Must run global alignment first.'); end
    % (A) IF ALIGNING TO A SCAN INDEX AND THEN SHIFTING GLOBALLY FOR ALL SCANS BY A CONSTANT
    if alignType == "global shift via scan"
        % - Extracting the scan to make the alignment to
        ebindx                      = scan_indxs(1);
        % - Eb alignment over a single defined scan
        [~, eb_shift, ~]            = AlignEF(dataStr.(dField)(:,:,ebindx), dataStr.eb(:,1,ebindx), eWin(1), dEWin, dESmooth, feat);
        % -- Applying the energy shift to the data globally
        dataStr.eb                  = dataStr.eb - eb_shift;
        % -- Appending all the energy shifts into a cell
        dataStr.eb_shifts{end+1}    = [ebindx, eb_shift];
    % (B) IF ALIGNING ONLY CERTAIN SCAN INDICES
    elseif alignType == "scans"
        ebindx      = 0;
        i_eWin    	= [];
        for i = scan_indxs
            ebindx = ebindx + 1;
            % - For a single value of eWin to be used
            if isempty(eWin) || length(eWin) == 1; i_eWin = eWin;
            % - For an eWin value at each scan between eWin(1) -> eWin(2)
            else; i_eWin = eWin(ebindx);
            end
            % - Eb alignment over defined scans
            [~, eb_shift, ~] = AlignEF(dataStr.(dField)(:,:,i), dataStr.eb(:,1,i), i_eWin, dEWin, dESmooth, feat);
            % -- Applying the energy shift to the data
            dataStr.eb(:,:,i) = dataStr.eb(:,:,i) - eb_shift;
            % -- Appending all the energy shifts into a matrix
            eb_shifts(i) = eb_shift;
        end
        dataStr.eb_shifts{end+1} = [scan_indxs, eb_shifts];
    end
end 
% - Appending new data to the data object
if ~isfield(dataStr.meta, 'align_args'); dataStr.meta.align_args = {}; end
dataStr.meta.align_args{end+1} = align_args;

%% - 3 - Remapping matrices into 2D or 3D consistent forms
% - 3.1 Remapping to 2D
% - Remapping the dataStr.eb domain to be 2D
if size(dataStr.eb,1) == 1 || size(dataStr.eb,2) == 1; dataStr.eb=repmat(dataStr.eb,[1, size(dataStr.(dField),2)]); end
% - 3.2 Remapping to 3D
% - Remapping the dataStr.eb domain to be 3D
if size(dataStr.eb, 3) == 1; dataStr.eb = repmat(dataStr.eb,[1, 1, size(dataStr.(dField),3)]); end
% - 3.3 Remap the tht matrix
if ~isfield(dataStr, 'tht')
    dataStr.tht = dataStr.raw_tht;
    % - 3.1 Remapping to 2D
    % -- Remapping the dataStr.tht / K domain to be 2D
    if size(dataStr.tht,1) == 1 || size(dataStr.tht,2) == 1; dataStr.tht=repmat(dataStr.tht,[size(dataStr.(dField),1), 1]); end
    % - 3.2 Remapping to 3D
    % -- Remapping the dataStr.tht / K domain to be 3D
    if size(dataStr.tht, 3) == 1; dataStr.tht = repmat(dataStr.tht,[1, 1, size(dataStr.(dField),3)]); end
elseif size(dataStr.tht, 1) == 1
    dataStr.tht = dataStr.raw_tht;
    % - 3.1 Remapping to 2D
    % -- Remapping the dataStr.tht / K domain to be 2D
    if size(dataStr.tht,1) == 1 || size(dataStr.tht,2) == 1; dataStr.tht=repmat(dataStr.tht,[size(dataStr.(dField),1), 1]); end
    % - 3.2 Remapping to 3D
    % -- Remapping the dataStr.tht / K domain to be 3D
    if size(dataStr.tht, 3) == 1; dataStr.tht = repmat(dataStr.tht,[1, 1, size(dataStr.(dField),3)]); end
end

end