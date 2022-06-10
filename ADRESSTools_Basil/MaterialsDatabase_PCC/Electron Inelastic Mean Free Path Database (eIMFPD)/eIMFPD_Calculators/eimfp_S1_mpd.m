function imfp = eimfp_S1_mpd(ke_dat, material)
% imfp = eimfp_S1_mpd(ke_dat, material)
%   Function that determines the electron inelastic mean free path (IMFP) 
%   based on the S1 equation described by M. P. Seah [1]. The IMFP here
%   additionally depends on the value of Z and this formalism is compatible
%   for elemental or binary materials. In this function, you can define the 
%   material as a string and it will look it up the relevant parameters in 
%   the Material Properties Database (MPD) ('MPD_PCC.mat') and determine 
%   the imfp using the S2 equation.
%
%   REFERENCE:
%   [1] M. P. Seah, “An accurate and simple universal curve for the 
%       energy-dependent electron inelastic mean free path,” 
%       Surf. Interface Anal., vol. 44, no. 4, pp. 497–503, 2012, 
%       doi: 10.1002/sia.4816.
%
%   IN:
%   -   ke_dat:     N×1 vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%   -   material:  	string of the material whose imfp is to be determined; e.g. "Si", "SiO2", "InAs", "Al2O3"...
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]

%% Default parameters (Parameters for Silicon)
if nargin < 2; material = "Si";  end
if isempty(material); material = "Si"; end
%% - 1 - Extracting the material parameters from the materials database
% - Extracting the material properties
material_props = get_mpd_props(material);
% - Extracting the material properties required for the S1 formalism
rho     = material_props.DENSITY;
Nv      = material_props.ELECT_VALENCY;
M       = material_props.ATOM_MASS;
Egap    = material_props.ELE_BGAP;
Z       = material_props.ATOM_ZNUM;
%% - 2 - Determination of the IMFP via S1 formalism
imfp = eimfp_S1(ke_dat, rho, M, Egap, Z);   % extract imfp in Angstroms
%% Ensuring the imfp is a column vector
if size(imfp, 2) >1; imfp = imfp'; end
end