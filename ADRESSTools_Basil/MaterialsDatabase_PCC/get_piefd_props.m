function piefd_props = get_piefd_props(element)
%  piefd_props = get_piefd_props(element)
%   This is a function that extracts all of the parameters from the
%   Photoionisation Energy and Fluorescence Database (PIEFD) for the element
%   defined here by the user. This database contains all the binding
%   energies, emission / fluorescence lines, emission / fluorescence
%   relative intensities and fluorescence yield.
%
%   IN:
%   -   element:            char/string of the element; e.g. "H", "He", "Si", "In"...
%
%   OUT:
%   -   piefd_props:        MATLAB table which contains all of the PIEFD parameters.

%% - 1 - Extracting the material parameters from the materials properties database
% -- Loading in the materials database
PIEFD_PCC         = load('PIEFD_PCC.mat'); PIEFD_PCC = PIEFD_PCC.PIEFD_PCC;
% -- Making the input material be in a 'char' format
element        = char(element); 
% -- Find the row-number / index for the material of choice
atom_symbols    = PIEFD_PCC.ATOM_SYMB;
idx             = find(strcmp(atom_symbols(:), element));
% -- Extract the properties of the desired material
piefd_props  = PIEFD_PCC(idx,:);
% -- Printing output text
% fprintf("LOADED DATA: ID: %i, Name: %s, Symbol: %s", piefd_props.ID, piefd_props.ATOM_NAME{1}, piefd_props.ATOM_SYMB{1});
end