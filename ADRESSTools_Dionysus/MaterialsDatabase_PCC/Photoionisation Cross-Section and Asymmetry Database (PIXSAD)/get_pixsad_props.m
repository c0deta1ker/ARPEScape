function [hv, xsect, asymmetry] = get_pixsad_props(element, corelevel, hvi)
% [hv, xsect, asymmetry] = get_pixsad_props(element, corelevel, hvi)
%   This is a function that extracts all of the parameters from the
%   Photoionisation Cross-Section and Asymmetry Database (PIXSAD) for the 
%   element defined here by the user. This database contains all the
%   photoionisation cross-sections and asymmetry parameters for elements 1
%   - 103, for orbitals that are probed with soft x-ray radiation (20 -
%   1500 eV). 
%
%   IN:
%   -   element:    	char/string of the element; e.g. "H", "He", "Si", "In"...
%   -   corelevel:      char/string of the core-level to be probed; e.g. "1s", "2s", "2p", "4d"...
%   -   hvi:            photon energy value (or range [min, max]) over which the cross-section and asymmetries will be determined.
%
%   OUT:
%   -   hv:             1×N vector (or table) of the interpolated / nearest photon energies to the input based on the data.
%   -   xsect:          1×N vector (or table) of the photoionisation cross-sections.
%   -   asymmetry:      1×N vector (or table) of the asymmetries.
%   NOTE: if any of the outputs are NaN, this means you are outside of the
%   range of photon energies for the defined data.

%% Default parameters
if nargin < 2; corelevel = [];  end
if nargin < 3; hvi = [];  end
if isempty(corelevel); corelevel = []; end
if isempty(hvi); hvi = []; end
% -- Sorting the photon energy to ascending order
hvi = sort(hvi);

%% - 1 - Loading the PIXS properties from our MATLAB database
PIXSAD_PCC	= load('PIXSAD_PCC.mat'); PIXSAD_PCC = PIXSAD_PCC.PIXSAD_PCC;

%% - 2 - Find the index of the element defined by the user
% -- Find the element
if ~isempty(element)
    % -- Find the row-number / index for the element of choice (case-insensitive)
    element  	= char(element);
    ele_indx 	= find(strcmpi([PIXSAD_PCC.ATOM_SYMB], element));
    if isempty(ele_indx)
        msg = 'Element could not be identified. Only use atomic-symbols for elements 1 - 103; H, He, Li, Be..., No, Lr'; error(msg)
    end
% -- If no element is defined, return an error
else
    msg = 'Element could not be identified. Only use atomic-symbols for elements 1 - 103; H, He, Li, Be..., No, Lr'; 
    error(msg)
end

%% - 3 - Extracting the photoionisation data
%% - A - When no core-level is defined, output all the data
if isempty(corelevel) && isempty(hvi)
    data        = PIXSAD_PCC.DataTable{ele_indx};
    hv          = data(:,1);
    xsect       = data(:,2:2:end);
    asymmetry   = data(:,3:2:end);
elseif isempty(corelevel)
    data     	= PIXSAD_PCC.DataTable{ele_indx};
    hv          = data.hv;
    [~, lbindx]	= min(abs(hv - min(hvi(:))));
    [~, ubindx]	= min(abs(hv - max(hvi(:))));
    % -- For a single photon energy entry
    if lbindx == ubindx
        hv          = data(lbindx,1);
        xsect       = data(lbindx,2:2:end);
        asymmetry   = data(lbindx,3:2:end);
    % -- For a multiple photon energies
    else
        hv          = data(lbindx:ubindx,1);
        xsect       = data(lbindx:ubindx,2:2:end);
        asymmetry   = data(lbindx:ubindx,3:2:end);
    end
else
    %% - A - Find the index of the core-level defined by the user
    element_cls	= PIXSAD_PCC.FileNames{ele_indx};
    % -- Find the core-level for the element of choice
    if ~isempty(corelevel) && length(corelevel) == 2
        corelevel  	= lower(char(corelevel));
        cl_indx    	= find(contains(lower(element_cls), corelevel));
        if isempty(cl_indx)
            msg = sprintf("Core-level not found. For element %s, only the following core-levels have data; ", element); 
            for i = 1:length(element_cls); msg = msg + string(element_cls{i}(end-6:end-5)) + ", "; end
            error(msg)
        end
        % -- Extracting the data from the PIXS database
        Data = PIXSAD_PCC.DataInt{ele_indx}{cl_indx};
    % -- If no core-level are defined, make user input one
    else
        msg = sprintf("Core-level not found. For element %s, only the following core-levels have data; ", element); 
        for i = 1:length(element_cls); msg = msg + string(element_cls{i}(end-6:end-5)) + ", "; end
        error(msg)
    end
    %% - B - Extract the photoionisaion cross-section and asymmetry for the chosen photon energy
    % -- Printing output text for the element and core-level selected
%     fprintf("PIXS DATA: Symbol: %s, Core-Level: %s", element, corelevel);
    % --- If photon energy is empty, default to all values
    if isempty(hvi)
        hv          = Data(:,1);
        xsect       = Data(:,2);
        asymmetry   = Data(:,3);
    % --- For a range of photon energies
    elseif length(hvi) > 1
        for i = 1:length(hvi)
            [~, hv_indx]	= min(abs(Data(:,1) - hvi(i)));
            % -- Extracting the relevant parameters
            hv(i)          = Data(hv_indx,1);
            xsect(i)       = Data(hv_indx,2);
            asymmetry(i)   = Data(hv_indx,3);
        end
    % --- For a single photon energy
    else
        % -- Find the index for the photon enegy
        [~, hv_indx]	= min(abs(Data(:,1) - hvi));
        % -- Extracting the relevant parameters
        hv          = Data(hv_indx,1);
        xsect       = Data(hv_indx,2);
        asymmetry   = Data(hv_indx,3);
    end
end
end