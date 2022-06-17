function dataStr = align_energy(dataStr, align_args, dataStr_ref, plot_result)
% dataStr = align_energy(dataStr, align_args, dataStr_ref, plot_result)
%   This is a function that will align the 'dataStr' ARPES data-structure 
%   to the Fermi-level, VB or CB state. It operates using the
%   AlignEF() function conditions, but can also load in a reference data file
%   'dataStr_ref' that it will use as the energy reference.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   [EAlign, EF, Fail] = AlignEF(Data, ECorr [,eWin] [,dEWin] [,dESmooth] [,feat]);
%
%   IN:
%   -   dataStr:            data structure of the ARPES data.
%   -   align_args:         {1×6} cell of {Type,Scans,eWin, dEWin, dESmooth, feat}.
%                               Type:       "global"
%                                         	"align2ref"
%                                       	"align2ref-perchannel"
%                                       	"fit2ef" || "fit2ef-95%" (default 95% window size)
%                                       	"fit2ef-5%"
%                                           "fit2ef-10%"
%                                           "fit2ef-25%"
%                                           "fit2ef-50%"
%                                           "fit2ef-75%"
%                                           "fit2Au4f" || "fit2Au4f-95%" (default 95% window size)
%                                           "global shift via scan"
%                                         	"scans"
%                                        	"none"
%                               Scans:      for "global" or "ref" it is empty []. For "scans", it defines the scan indices to perform alignment on.
%                               eWin:       approximate EF position. Can be a single, constant value. If a vector, it linear interpolates estimates between [firstEF, endEF].
%                               dEWin:      can be empty []. Single, constant value that defines the width of the EF.
%                               dESmooth:   can be empty []. Single, constant value that defines Gaussian pre-smoothing energy width.
%                               feat: either 'edge' or 'peak'. defines the feature to align as a 'peak' instead of the default Fermi 'edge'.
%   -   dataStr_ref:     	reference data structure used for Eb alignment.
%   -   plot_result:        if 1, will plot figure of the alignment, otherwise it wont.
%
%   OUT:
%   dataStr - MATLAB data structure with new additional fields below;
%   -   .meta.align_args:  {1×6} cell of input arguments of alignment.
%	-   .(tht):         aligned 2D or 3D array of theta.
%	-   .(eb):          aligned 2D or 3D array of energy.
%	-   .(eb_shifts): 	cell array that contains all of the Eb shifts applied to this data.

%% Default parameters
if nargin < 4; plot_result = 0; end
if nargin < 3; dataStr_ref = dataStr; end
if nargin < 2; align_args = {"global", [], 0, [], [], 'edge'}; end
if isempty(align_args); align_args = {"global", [], 0, [], [], 'edge'}; end
if isempty(dataStr_ref); dataStr = dataStr_ref; end
if isempty(plot_result); plot_result = 0; end
% - Extracting the fields to be used with most recent processing
[~, yField, ~, dField]                  = find_data_fields(dataStr);
[xField_ref, yField_ref, ~, dField_ref] = find_data_fields(dataStr_ref);
%% Validity checks on the input parameters
if length(align_args) ~= 6
    error('error: align_args is not a {1x6} cell - make sure you defined the align arguments correctly.');
else
    align_args{1} = lower(string(align_args{1}));
    if isempty(align_args{2}) || ischar(align_args{2}) || isstring(align_args{2}); align_args{2} = []; end
    if isempty(align_args{3}) || ischar(align_args{3}) || isstring(align_args{3}); align_args{3} = []; end
    if isempty(align_args{4}) || ischar(align_args{4}) || isstring(align_args{4}); align_args{4} = []; end
    if isempty(align_args{5}) || ischar(align_args{5}) || isstring(align_args{5}); align_args{5} = []; end
    align_args{6} = char(lower(string(align_args{6})));
end

%% - 1 - Initialising the alignment parameters
% - Extracting alignment alignType parameters
alignType   = align_args{1};
scan_indxs	= align_args{2};
% - Extracting alignment parameters
eWin        = align_args{3};
dEWin       = align_args{4};
dESmooth    = align_args{5};
feat        = align_args{6};
% - Defining the scan indices over the energy alignment for-loop
if ~strcmpi(alignType,"scans") && ~strcmpi(alignType,"global shift via scan")
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
if isempty(eWin);           eWin = [];
elseif length(eWin) == 1;   eWin = ones(length(scan_indxs)) .* eWin;
elseif length(eWin) > 1;    eWin = linspace(eWin(1), eWin(2), length(scan_indxs));
end

%% - 2 - Executing binding energy alignment over all scans
% - Checking the align type is valid
if strcmpi(alignType,"global") || strcmpi(alignType,"global shift via scan") || strcmpi(alignType,"scans")...
        || strcmpi(alignType,"align2ref") || strcmpi(alignType,"align2ref-perchannel")...
        || strcmpi(alignType,"fit2ef") || strcmpi(alignType,"fit2ef-5%") || strcmpi(alignType,"fit2ef-10%") || strcmpi(alignType,"fit2ef-25%")...
        || strcmpi(alignType,"fit2ef-50%") || strcmpi(alignType,"fit2ef-75%") || strcmpi(alignType,"fit2ef-95%")...
        || strcmpi(alignType,"fit2Au4f") || strcmpi(alignType,"fit2Au4f-95%")...
        || strcmpi(alignType,"none")
else
    error('Error: Align type not recognised.');
end
%% EXECUTING GLOBAL ENERGY ALIGNMENTS
if strcmpi(alignType,"global")...
        || strcmpi(alignType,"align2ref") || strcmpi(alignType,"align2ref-perchannel")...
        || strcmpi(alignType,"fit2ef") || strcmpi(alignType,"fit2ef-5%") || strcmpi(alignType,"fit2ef-10%") || strcmpi(alignType,"fit2ef-25%") || strcmpi(alignType,"fit2ef-50%") || strcmpi(alignType,"fit2ef-75%") || strcmpi(alignType,"fit2ef-95%")
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
        if strcmpi(alignType,"global")
            % --- If the alignment was previously performed, apply again
            if string(yField) == "eb";          [~, eb_shift, ~] = AlignEF(dataStr.(dField)(:,:,i), dataStr.eb(:,1,i), i_eWin, dEWin, dESmooth, feat);
            % --- If no alignment has been performed, do it for the first time
            elseif string(yField) == "raw_eb";  [~, eb_shift, ~] = AlignEF(dataStr.(dField)(:,:,i), dataStr.raw_eb, i_eWin, dEWin, dESmooth, feat);
            end
        % -- (B) IF ALIGNING GLOBALLY TO A REFERENCE DATA-SET
        elseif strcmpi(alignType,"align2ref")
            % --- If the scan index exceeds the reference data, cap it to the end
            if i > size(dataStr_ref.(dField_ref), 3)
                [~, eb_shift, ~] = AlignEF(dataStr_ref.(dField_ref)(:,:,end), dataStr_ref.(yField_ref), i_eWin, dEWin, dESmooth, feat);
            else
                [~, eb_shift, ~] = AlignEF(dataStr_ref.(dField_ref)(:,:,i), dataStr_ref.(yField_ref), i_eWin, dEWin, dESmooth, feat);
            end
        % -- (C) IF ALIGNING GLOBALLY TO A REFERENCE DATA-SET, TO EACH INDIVIDUAL X-CHANNEL
        elseif strcmpi(alignType,"align2ref-perchannel")
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
        elseif strcmpi(alignType,"fit2ef") || strcmpi(alignType,"fit2ef-95%")...
                || strcmpi(alignType,"fit2ef-5%")...
                || strcmpi(alignType,"fit2ef-10%")...
                || strcmpi(alignType,"fit2ef-25%")...
                || strcmpi(alignType,"fit2ef-50%")...
                || strcmpi(alignType,"fit2ef-75%")
            % --- Initialising variables
            if strcmpi(alignType,"fit2ef") || strcmpi(alignType,"fit2ef-95%");	angleWin = mean(dataStr_ref.(xField_ref)(:)) + [-1,1].*0.95*abs(range(dataStr_ref.(xField_ref)(:)));
            elseif strcmpi(alignType,"fit2ef-25%");  angleWin = mean(dataStr_ref.(xField_ref)(:)) + [-1,1].*0.25*abs(range(dataStr_ref.(xField_ref)(:)));
            elseif strcmpi(alignType,"fit2ef-50%");  angleWin = mean(dataStr_ref.(xField_ref)(:)) + [-1,1].*0.50*abs(range(dataStr_ref.(xField_ref)(:)));
            elseif strcmpi(alignType,"fit2ef-75%");  angleWin = mean(dataStr_ref.(xField_ref)(:)) + [-1,1].*0.75*abs(range(dataStr_ref.(xField_ref)(:)));
            elseif strcmpi(alignType,"fit2ef-5%");   angleWin = mean(dataStr_ref.(xField_ref)(:)) + [-1,1].*0.05*abs(range(dataStr_ref.(xField_ref)(:)));
            elseif strcmpi(alignType,"fit2ef-10%");  angleWin = mean(dataStr_ref.(xField_ref)(:)) + [-1,1].*0.10*abs(range(dataStr_ref.(xField_ref)(:)));
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
        % -- (E) IF ALIGNING TO THE Au4f CORE-LEVEL
        elseif strcmpi(alignType,"fit2Au4f") || strcmpi(alignType,"fit2Au4f-95%")
            Au4f72 = -84.00;
            [pesStr, ~] = pes_extract_data(arpesStr, {angleWin,1}, 0);
            [xVal, ~]   = find_peak_loc(pesStr.xdat, pesStr.ydat, Au4f72 + [-0.5, 0.5], "sGLA", plot_result);
            eb_shift    = xVal - Au4f72;
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
if strcmpi(alignType,"global shift via scan") || strcmpi(alignType,"scans")
    if ~isfield(dataStr, 'eb'); error('Must run global alignment first.'); end
    % (A) IF ALIGNING TO A SCAN INDEX AND THEN SHIFTING GLOBALLY FOR ALL SCANS BY A CONSTANT
    if strcmpi(alignType,"global shift via scan")
        % - Extracting the scan to make the alignment to
        ebindx                      = scan_indxs(1);
        % - Eb alignment over a single defined scan
        [~, eb_shift, ~]            = AlignEF(dataStr.(dField)(:,:,ebindx), dataStr.eb(:,1,ebindx), eWin(1), dEWin, dESmooth, feat);
        % -- Applying the energy shift to the data globally
        dataStr.eb                  = dataStr.eb - eb_shift;
        % -- Appending all the energy shifts into a cell
        dataStr.eb_shifts{end+1}    = [ebindx, eb_shift];
    % (B) IF ALIGNING ONLY CERTAIN SCAN INDICES
    elseif strcmpi(alignType,"scans")
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