function [dataStr, dataStr_ref] = align_energy(dataStr, align_args, dataStr_ref, plot_result)
% [dataStr, dataStr_ref] = align_energy(dataStr, align_args, dataStr_ref, plot_result)
%   - #1 - FIRST STEP IN PROCESSING ARPES DATA
%   This is a function that will align the 'dataStr' ARPES data-structure 
%   to the Fermi-level, a peak, or reference VB or CB state. The energy
%   alignment is done using various approaches that the user can select from,
%   either using the 'dataStr' or a reference 'dataStr_ref' data file. 
%       (1) 'align2ref'  : aligns to some peak or edge reference using MATools AlignEF() function.
%       (2) 'fit2ef'    : fits directly to the FDD edge of a reference. 
%       (3) 'fit2peak'  : fits directly to some peak feature.
%       (4) 'gsv1s'     : global shift via one scan aligned to some peak / edge reference using MATools AlignEF() function.
%
%   IN:
%   -   dataStr:            data structure of the ARPES data.
%   -   align_args:         {1×7} cell of {Type,Scans,eWin, dEWin, dESmooth, feat}.
%                               Type:       "align2ref" / "align2ef" / "alignef" / "global" (same as "AlignEF")
%                                       	"fit2ef" / "fit2FDD" 
%                                           "fit2FDDG" or "fit2FDDGpL" or "fit2FDDGsL"
%                                           "fit2peak"
%                                           "gsv1s" or "global shift via scan"
%                                        	"none"
%                               Scans:      can be empty [] to align energy to all scans. Otherwise, you can define a list of the scan indices to align energy.
%                               eWin:       approximate EF position. Can be a single, constant value. If a vector, it linear interpolates estimates between [firstEF, endEF].
%                               dEWin:      can be empty []. Single, constant value that defines the width of the EF.
%                               dESmooth:   can be empty []. Single, constant value that defines Gaussian pre-smoothing energy width.
%                               feat:       either 'edge' or 'peak'. defines the feature to align as a 'peak' instead of the default Fermi 'edge'.
%                               kWin:       a constant value that defines the angle/momentum window to use as a % of the full window. Default: 0.95.
%   -   dataStr_ref:     	reference data structure used for Eb alignment.
%   -   plot_result:        if 1, will plot figure of the alignment, otherwise it wont.
%
%   OUT:
%   dataStr - MATLAB data structure with new additional fields below;
%   -   .meta.align_args:   {1×7} cell of input arguments of alignment.
%	-   .(eb_shifts): 	    cell array that contains all of the Eb shifts applied to this data.
%	-   .(eb_ref_file):     cell array that contains the name of Eb reference files used.
%	-   .(eb):              aligned 2D or 3D array of energy.
%	-   .(tht):             aligned 2D or 3D array of theta.
%   dataStr_ref - MATLAB reference data structure with new additional fields below;
%   -   .meta.align_args:   {1×7} cell of input arguments of alignment.
%	-   .(eb_shifts): 	    cell array that contains all of the Eb shifts applied to this data.
%	-   .(raw_eb):          aligned 2D array of energy.

%% Default parameters
if nargin < 4; plot_result = 0; end
if nargin < 3; dataStr_ref = dataStr; end
if nargin < 2; align_args = {"align2ref", [], 0, [], [], "edge", 0.95}; end
if isempty(align_args); align_args = {"align2ref", [], 0, [], [], "edge", 0.95}; end
if isempty(dataStr_ref); dataStr_ref = dataStr; end
if isempty(plot_result); plot_result = 0; end
% - Extracting the fields to be used with most recent processing
[xField, yField, ~, dField]             = find_data_fields(dataStr);
[xField_ref, yField_ref, ~, dField_ref] = find_data_fields(dataStr_ref);
% - Constants
eb_ref_file     = dataStr_ref.FileName;

%% - 1 - Validity check on the input parameters
% -- Checking the length of input arguments
if length(align_args) == 1; align_args = {align_args{1}, [], [], [], [], [], []};
elseif length(align_args) == 2; align_args = {align_args{1}, align_args{2}, [], [], [], [], []};
elseif length(align_args) == 3; align_args = {align_args{1}, align_args{2}, align_args{3}, [], [], [], []};
elseif length(align_args) == 4; align_args = {align_args{1}, align_args{2}, align_args{3}, align_args{4}, [], [], []};
elseif length(align_args) == 5; align_args = {align_args{1}, align_args{2}, align_args{3}, align_args{4}, align_args{5}, [], []};
elseif length(align_args) == 6; align_args = {align_args{1}, align_args{2}, align_args{3}, align_args{4}, align_args{5}, align_args{6}, []};
elseif length(align_args) > 7; error('align_args is not a {1x7} cell. Make sure you defined the align arguments correctly. Index must not exceed 7.');
end
%% - 2 - Initialising the energy alignment variables
% -- Default values if undefined
alignType   = align_args{1}; if isempty(alignType);     alignType = "align2ef"; end
sIndxs      = align_args{2}; if isempty(sIndxs);        sIndxs    = []; end
eWin        = align_args{3}; if isempty(eWin);          eWin      = 0.0; end
dEWin       = align_args{4}; if isempty(dEWin);         dEWin     = []; end
dESmooth    = align_args{5}; if isempty(dESmooth);      dESmooth  = []; end
feat        = align_args{6}; if isempty(feat);          feat      = "edge"; end
kWin        = align_args{7}; if isempty(kWin);          kWin      = 0.95; end
% -- Making all string lowercase
alignType   = lower(string(alignType));
feat        = lower(string(feat));
% -- Checking the align type is valid
if strcmpi(alignType,"align2ref") || strcmpi(alignType,"align2ef") ||strcmpi(alignType,"alignef") || strcmpi(alignType,"global")...
        || strcmpi(alignType,"fit2ef") || strcmpi(alignType,"fit2FDD") || strcmpi(alignType,"fit2FDDG") || strcmpi(alignType,"fit2FDDGpL") || strcmpi(alignType,"fit2FDDGsL")...
        || strcmpi(alignType,"fit2peak")...
        || strcmpi(alignType,"gsv1s") || strcmpi(alignType,"global shift via scan")...
        || strcmpi(alignType,"none")...
        || contains(alignType,"fit2ef") || contains(alignType,"fit2peak")
else; error('Error: Align Type not recognised.');
end
% -- Defining the scan indices over the energy alignment for-loop
if isempty(sIndxs) || ~strcmpi(alignType,"gsv1s") && ~strcmpi(alignType,"global shift via scan") && ~strcmpi(alignType,"scans")
    sIndxs = 1:size(dataStr.(dField), 3);
else
    minVal = 1; maxVal = size(dataStr.(dField), 3);
    % - Checking max/min are not exceeded
    sIndxs(sIndxs < minVal) = minVal;
    sIndxs(sIndxs > maxVal) = maxVal;
    % - Checking that there are not duplicate numbers
    sIndxs = unique(sIndxs);
end
% - Defining the values of eWin and interpolating if necessary
if isempty(eWin);           eWin = [];
elseif length(eWin) == 1;   eWin = ones(length(sIndxs)) .* eWin;
elseif length(eWin) > 1;    eWin = linspace(eWin(1), eWin(2), length(sIndxs));
end

%% - 3.1 - EXECUTING GLOBAL ENERGY ALIGNMENTS (filing through all scan indices)
if strcmpi(alignType,"align2ref") || strcmpi(alignType,"align2ef") ||strcmpi(alignType,"alignef") || strcmpi(alignType,"global")...
        || strcmpi(alignType,"fit2ef") || strcmpi(alignType,"fit2FDD") || strcmpi(alignType,"fit2FDDG") || strcmpi(alignType,"fit2FDDGpL") || strcmpi(alignType,"fit2FDDGsL")...
        || strcmpi(alignType,"fit2peak")...
        || contains(alignType,"fit2ef") || contains(alignType,"fit2peak")
    % - Initialising variables
    eb_shift_mat       = [];
    aligned_energy	= [];
    eWin_i          = [];
    % - Iterating through all scan indices
    for i = sIndxs
        % -- For a single value of eWin to be used
        if isempty(eWin) || length(eWin) == 1; eWin_i = eWin;
        % -- For an eWin value at each scan between eWin(1) -> eWin(2)
        else; eWin_i = eWin(i);
        end

        % -- (A) IF ALIGNING GLOBALLY TO A REFERENCE DATA-SET
        if strcmpi(alignType,"global") || strcmpi(alignType,"align2ref") || strcmpi(alignType,"align2ef") || strcmpi(alignType,"alignef")
            % --- If the alignment was previously performed, apply again
            if string(yField_ref) == "eb"
                % --- If the scan index exceeds the reference data, cap it to the end
                if i > size(dataStr_ref.(dField_ref), 3);   [~, eb_shift, ~] = AlignEF(dataStr_ref.(dField_ref)(:,:,end), dataStr_ref.eb(:,1,end), eWin_i, dEWin, dESmooth, feat);
                else;                                       [~, eb_shift, ~] = AlignEF(dataStr_ref.(dField_ref)(:,:,i), dataStr_ref.eb(:,1,i), eWin_i, dEWin, dESmooth, feat);
                end
            % --- If no alignment has been performed, do it for the first time
            elseif string(yField_ref) == "raw_eb"
                % --- If the scan index exceeds the reference data, cap it to the end
                if i > size(dataStr_ref.(dField_ref), 3);   [~, eb_shift, ~] = AlignEF(dataStr_ref.(dField_ref)(:,:,end), dataStr_ref.raw_eb, eWin_i, dEWin, dESmooth, feat);
                else;                                       [~, eb_shift, ~] = AlignEF(dataStr_ref.(dField_ref)(:,:,i), dataStr_ref.raw_eb, eWin_i, dEWin, dESmooth, feat);
                end
            end

        % -- (B) IF FITTING TO THE FERMI-EDGE OF A REFERENCE DATA-SET
        elseif strcmpi(alignType,"fit2ef") || strcmpi(alignType,"fit2FDDGsL") || contains(alignType,"fit2ef")
            % --- Initialising variables
            angleWin = mean(dataStr_ref.(xField_ref)(:)) + [-1,1].*kWin*abs(range(dataStr_ref.(xField_ref)(:)));
            % ---- Executing energy alignment by fitting to FDD
            if i > size(dataStr_ref.raw_data, 3);       [~,eb_shift] = fdd2fit_solver(dataStr_ref.(xField_ref), dataStr_ref.(yField_ref), dataStr_ref.(dField_ref)(:,:,end), angleWin, eWin_i+dEWin.*[-1,1], "FDDGsL", [], [], plot_result);
            else;                                       [~,eb_shift] = fdd2fit_solver(dataStr_ref.(xField_ref), dataStr_ref.(yField_ref), dataStr_ref.(dField_ref)(:,:,i), angleWin, eWin_i+dEWin.*[-1,1], "FDDGsL", [], [], plot_result);
            end
        elseif strcmpi(alignType,"fit2FDD")
            % --- Initialising variables
            angleWin = mean(dataStr_ref.(xField_ref)(:)) + [-1,1].*kWin*abs(range(dataStr_ref.(xField_ref)(:)));
            % ---- Executing energy alignment by fitting to FDD
            if i > size(dataStr_ref.raw_data, 3);       [~,eb_shift] = fdd2fit_solver(dataStr_ref.(xField_ref), dataStr_ref.(yField_ref), dataStr_ref.(dField_ref)(:,:,end), angleWin, eWin_i+dEWin.*[-1,1], "FDD", [], [], plot_result);
            else;                                       [~,eb_shift] = fdd2fit_solver(dataStr_ref.(xField_ref), dataStr_ref.(yField_ref), dataStr_ref.(dField_ref)(:,:,i), angleWin, eWin_i+dEWin.*[-1,1], "FDD", [], [], plot_result);
            end
        elseif strcmpi(alignType,"fit2FDDG")
            % --- Initialising variables
            angleWin = mean(dataStr_ref.(xField_ref)(:)) + [-1,1].*kWin*abs(range(dataStr_ref.(xField_ref)(:)));
            % ---- Executing energy alignment by fitting to FDD
            if i > size(dataStr_ref.raw_data, 3);       [~,eb_shift] = fdd2fit_solver(dataStr_ref.(xField_ref), dataStr_ref.(yField_ref), dataStr_ref.(dField_ref)(:,:,end), angleWin, eWin_i+dEWin.*[-1,1], "FDDG", [], [], plot_result);
            else;                                       [~,eb_shift] = fdd2fit_solver(dataStr_ref.(xField_ref), dataStr_ref.(yField_ref), dataStr_ref.(dField_ref)(:,:,i), angleWin, eWin_i+dEWin.*[-1,1], "FDDG", [], [], plot_result);
            end
        elseif strcmpi(alignType,"fit2FDDGpL")
            % --- Initialising variables
            angleWin = mean(dataStr_ref.(xField_ref)(:)) + [-1,1].*kWin*abs(range(dataStr_ref.(xField_ref)(:)));
            % ---- Executing energy alignment by fitting to FDD
            if i > size(dataStr_ref.raw_data, 3);       [~,eb_shift] = fdd2fit_solver(dataStr_ref.(xField_ref), dataStr_ref.(yField_ref), dataStr_ref.(dField_ref)(:,:,end), angleWin, eWin_i+dEWin.*[-1,1], "FDDGpL", [], [], plot_result);
            else;                                       [~,eb_shift] = fdd2fit_solver(dataStr_ref.(xField_ref), dataStr_ref.(yField_ref), dataStr_ref.(dField_ref)(:,:,i), angleWin, eWin_i+dEWin.*[-1,1], "FDDGpL", [], [], plot_result);
            end

        % -- (C) IF ALIGNING TO A PEAK FEATURE OF A REFERENCE DATA-SET
        elseif strcmpi(alignType,"fit2peak") || contains(alignType,"fit2peak")
            % --- Initialising variables
            angleWin = mean(dataStr_ref.(xField_ref)(:)) + [-1,1].*kWin*abs(range(dataStr_ref.(xField_ref)(:)));
            % --- Extracting the peak position
            [pesStr, ~] = pes_extract_data(dataStr_ref, {angleWin,1}, 0);
            if isempty(eWin);       [xVal, ~]   = find_peak_loc(pesStr.xdat, pesStr.ydat, [], "spline", plot_result);
            elseif isempty(dEWin);  [xVal, ~]   = find_peak_loc(pesStr.xdat, pesStr.ydat, [], "spline", plot_result);
            else;                   [xVal, ~]   = find_peak_loc(pesStr.xdat, pesStr.ydat, eWin + dEWin.*[-1, 1], "spline", plot_result);
            end
            eb_shift    = xVal - eWin;
        end

        % - Appending all the aligned data vectors into a matrix
        % -- Actual Data
        if isfield(dataStr, 'eb');  aligned_energy(:,:,i) = dataStr.eb(:,:,i) - eb_shift;
        else;                       aligned_energy(:,1,i) = dataStr.raw_eb - eb_shift;
        end
        % -- Reference Data (ensuring it keeps the same size)
        if i <= size(dataStr_ref.(dField_ref), 3)
            if isfield(dataStr_ref, 'eb');  aligned_energy_ref(:,:,i) = dataStr_ref.eb(:,:,i) - eb_shift;
            else;                           aligned_energy_ref(:,1,i) = dataStr_ref.raw_eb - eb_shift;
            end
        end
        % - Appending all the energy shifts into a matrix
        eb_shift_mat(i,:) = [i, eb_shift];
    end

    % Assigning final variables to the data structure
    EB      = aligned_energy;
    EB_REF  = aligned_energy_ref;
end

%% - 3.2 - EXECUTING SCAN SELECTED ENERGY ALIGNMENTS (filing through selected scan indices)
if strcmpi(alignType,"global shift via scan") || strcmpi(alignType,"scans")
    if ~isfield(dataStr, 'eb'); error('Must run global alignment first.'); end

    % (A) IF ALIGNING TO A SCAN INDEX AND THEN SHIFTING GLOBALLY FOR ALL SCANS BY A CONSTANT
    if strcmpi(alignType,"global shift via scan") || strcmpi(alignType,"gsv1s")
        % - Extracting the scan to make the alignment to
        ebindx              = sIndxs(1);
        % - Eb alignment over a single defined scan
        [~, eb_shift, ~]    = AlignEF(dataStr_ref.(dField_ref)(:,:,ebindx), dataStr_ref.eb(:,1,ebindx), eWin(1), dEWin, dESmooth, feat);
        % -- Applying the energy shift to the data globally
        EB                  = dataStr.eb - eb_shift;
        EB_REF              = dataStr_ref.eb - eb_shift;
        % -- Appending all the energy shifts into a cell
        eb_shift_mat        = [ebindx, eb_shift];
        
    % (B) IF ALIGNING ONLY CERTAIN SCAN INDICES
    elseif strcmpi(alignType,"scans")
        ebindx      = 0;
        eWin_i    	= [];
        eb_shift_mat = [];
        EB      = dataStr.eb;
        EB_REF  = dataStr_ref.eb;
        for i = sIndxs
            ebindx = ebindx + 1;
            % - For a single value of eWin to be used
            if isempty(eWin) || length(eWin) == 1; eWin_i = eWin;
            % - For an eWin value at each scan between eWin(1) -> eWin(2)
            else; eWin_i = eWin(ebindx);
            end
            % - Eb alignment over defined scans
            [~, eb_shift, ~]    = AlignEF(dataStr_ref.(dField_ref)(:,:,i), dataStr_ref.eb(:,1,i), eWin_i, dEWin, dESmooth, feat);
            % -- Applying the energy shift to the data
            EB(:,:,i)                  = dataStr.eb(:,:,i) - eb_shift;
            EB_REF(:,:,i)              = dataStr_ref.eb(:,:,i) - eb_shift;
            % -- Appending all the energy shifts into a matrix
            eb_shift_mat(end+1) = [i, eb_shift];
        end
    end
end 

%% - 3 - Remapping matrices into 2D or 3D consistent forms
%% ARPES Data
% - 3.1 Remapping to 2D
% - Remapping the EB domain to be 2D
if size(EB,1) == 1 || size(EB,2) == 1; EB=repmat(EB,[1, size(dataStr.(dField),2)]); end
% - 3.2 Remapping to 3D
% - Remapping the EB domain to be 3D
if size(EB, 3) == 1; EB = repmat(EB,[1, 1, size(dataStr.(dField),3)]); end
% - 3.3 Remap the tht matrix
if ~isfield(dataStr, 'tht')
    THT = dataStr.raw_tht;
    % - 3.1 Remapping to 2D
    % -- Remapping the THT / K domain to be 2D
    if size(THT,1) == 1 || size(THT,2) == 1; THT=repmat(THT,[size(dataStr.(dField),1), 1]); end
    % - 3.2 Remapping to 3D
    % -- Remapping the THT / K domain to be 3D
    if size(THT, 3) == 1; THT = repmat(THT,[1, 1, size(dataStr.(dField),3)]); end
elseif size(dataStr.tht, 1) == 1
    THT = dataStr.raw_tht;
    % - 3.1 Remapping to 2D
    % -- Remapping the THT / K domain to be 2D
    if size(THT,1) == 1 || size(THT,2) == 1; THT=repmat(THT,[size(dataStr.(dField),1), 1]); end
    % - 3.2 Remapping to 3D
    % -- Remapping the THT / K domain to be 3D
    if size(THT, 3) == 1; THT = repmat(THT,[1, 1, size(dataStr.(dField),3)]); end
else; THT = dataStr.tht;
end

%% Reference Data
% - 3.1 Remapping to 2D
% - Remapping the EB_REF domain to be 2D
if size(EB_REF,1) == 1 || size(EB_REF,2) == 1; EB_REF=repmat(EB_REF,[1, size(dataStr_ref.(dField),2)]); end
% - 3.2 Remapping to 3D
% - Remapping the EB_REF domain to be 3D
if size(EB_REF, 3) == 1; EB_REF = repmat(EB_REF,[1, 1, size(dataStr_ref.(dField),3)]); end
% - 3.3 Remap the THT_REF matrix
if ~isfield(dataStr_ref, 'tht')
    THT_REF = dataStr_ref.raw_tht;
    % - 3.1 Remapping to 2D
    % -- Remapping the THT_REF / K domain to be 2D
    if size(THT_REF,1) == 1 || size(THT_REF,2) == 1; THT_REF=repmat(THT_REF,[size(dataStr_ref.(dField),1), 1]); end
    % - 3.2 Remapping to 3D
    % -- Remapping the THT_REF / K domain to be 3D
    if size(THT_REF, 3) == 1; THT_REF = repmat(THT_REF,[1, 1, size(dataStr_ref.(dField),3)]); end
elseif size(dataStr_ref.tht, 1) == 1
    THT_REF = dataStr_ref.raw_tht;
    % - 3.1 Remapping to 2D
    % -- Remapping the THT_REF / K domain to be 2D
    if size(THT_REF,1) == 1 || size(THT_REF,2) == 1; THT_REF=repmat(THT_REF,[size(dataStr_ref.(dField),1), 1]); end
    % - 3.2 Remapping to 3D
    % -- Remapping the THT_REF / K domain to be 3D
    if size(THT_REF, 3) == 1; THT_REF = repmat(THT_REF,[1, 1, size(dataStr_ref.(dField),3)]); end
else; THT_REF = dataStr_ref.tht;
end

%% 4 - Appending the new data to the MATLAB structure
%% ARPES Data
if ~isfield(dataStr.meta, 'align_args'); dataStr.meta.align_args = {}; end
dataStr.meta.align_args{end+1} = align_args;
if ~isfield(dataStr, 'eb_shifts'); dataStr.eb_shifts = {}; end
dataStr.eb_shifts{end+1}    = eb_shift_mat;
if ~isfield(dataStr, 'eb_ref_file'); dataStr.eb_ref_file = {}; end
dataStr.eb_ref_file{end+1} = eb_ref_file;
dataStr.eb  = EB;
dataStr.tht = THT;
%% Reference Data
if ~isfield(dataStr_ref.meta, 'align_args'); dataStr_ref.meta.align_args = {}; end
dataStr_ref.meta.align_args{end+1} = align_args;
if ~isfield(dataStr_ref, 'eb_shifts'); dataStr_ref.eb_shifts = {}; end
dataStr_ref.eb_shifts{end+1}    = eb_shift_mat;
dataStr_ref.eb  = EB_REF;
dataStr_ref.tht = THT_REF;

%% -- For Debugging
if plot_result == 1
    % - Extracting the fields to be used with most recent processing
    [xField, yField, ~, dField]             = find_data_fields(dataStr);
    [xField_ref, yField_ref, ~, dField_ref] = find_data_fields(dataStr_ref);
    pp = plot_props();
    fig = figure(); 
    fig.Position(3) = 2*pp.fig5x4(1); 
    fig.Position(4) = pp.fig5x4(2);
    subplot(121); hold on;
    ImData(dataStr_ref.(xField_ref)(:,:,sIndxs(1)), dataStr_ref.(yField_ref)(:,:,sIndxs(1)), dataStr_ref.(dField_ref)(:,:,sIndxs(1)));
    yline(0, 'color', 'c', 'linewidth', 1.5, 'linestyle', '-');
    yline(-1*eb_shift_mat(2), 'color', 'c', 'linewidth', 0.5, 'linestyle', '--'); 
    img_props(0, 'Color', [1 1 1]); colormap(hot); 
    axis([min(dataStr_ref.(xField_ref)(:)), max(dataStr_ref.(xField_ref)(:)), min(dataStr_ref.(yField_ref)(:)), max(dataStr_ref.(yField_ref)(:))]);
    colorbar; title('align_energy(): Reference', 'Interpreter', 'none');
    text(0.05, 0.25, "$$ \Delta E = $$ " + string(round(eb_shift_mat(2),3)) + " eV", 'interpreter', 'latex', 'fontsize', 12, 'color', 'w', 'Units','normalized');
    subplot(122); hold on;
    ImData(dataStr.(xField)(:,:,sIndxs(1)), dataStr.(yField)(:,:,sIndxs(1)), dataStr.(dField)(:,:,sIndxs(1)));
    yline(0, 'color', 'c', 'linewidth', 1.5, 'linestyle', '-');
    yline(-1*eb_shift_mat(2), 'color', 'c', 'linewidth', 0.5, 'linestyle', '--'); 
    img_props(0, 'Color', [1 1 1]); colormap(hot); 
    axis([min(dataStr.(xField)(:)), max(dataStr.(xField)(:)), min(dataStr.(yField)(:)), max(dataStr.(yField)(:))]);
    colorbar; title('align_energy(): Data', 'Interpreter', 'none');
    text(0.05, 0.25, "$$ \Delta E = $$ " + string(round(eb_shift_mat(2),3)) + " eV", 'interpreter', 'latex', 'fontsize', 12, 'color', 'w', 'Units','normalized');
end
end