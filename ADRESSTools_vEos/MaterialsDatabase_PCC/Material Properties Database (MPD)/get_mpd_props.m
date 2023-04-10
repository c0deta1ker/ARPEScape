function material_props = get_mpd_props(material)
% material_props = get_mpd_props(material)
%   This is a function that extracts all of the parameters from the
%   Materials Properties Database (MPD) for the material defined here by
%   the user.
%
%   IN:
%   -   material:           char/string of the material; e.g. "Si", "SiO2", "Al2O3"...
%
%   OUT:
%   -   material_props:     MATLAB table which contains all of the material parameters.

%% Default parameters
if nargin < 1; material = []; end
if isempty(material); material = []; end

%% - 1 - Extracting the material parameters from the materials properties database
% -- Loading in the materials database
MPD_PCC         = load('MPD_PCC.mat'); MPD_PCC = MPD_PCC.MPD_PCC;
% -- Making the input material be in a 'char' format
material        = char(material); 
% -- Find the row-number / index for the material of choice
atom_symbols    = MPD_PCC.ATOM_SYMB;
idx             = find(strcmpi(atom_symbols(:), material));
% -- Extract the properties of the desired material
if ~isempty(idx); material_props  = MPD_PCC(idx,:);
% -- If no material is identified, return an error
else; msg = 'Material could not be identified or does not exist in database.'; material_props = [];
end

end