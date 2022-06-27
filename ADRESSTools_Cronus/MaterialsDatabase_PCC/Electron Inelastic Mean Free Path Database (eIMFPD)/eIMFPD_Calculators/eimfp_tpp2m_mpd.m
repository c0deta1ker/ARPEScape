function imfp = eimfp_tpp2m_mpd(ke_dat, material)
% imfp = eimfp_tpp2m_mpd(ke_dat, material)
%   Function that determines the electron inelastic mean free path (IMFP) 
%   based on the TPP-2M energy compensation factor, which is a more
%   accurate method than the standard universal curve determination. The
%   input parameters of the materials however are now required. In this
%   function, you can define the material as a string and it will look it
%   up the relevant parameters in the Material Properties Database (MPD)
%   ('MPD_PCC.mat') and determine the imfp using the TPP-2M equation.
%
%   REFERENCE:
%       Tanuma, S., Powell, C. J., & Penn, D. R. (1987). 
%       Proposed formula for electron inelastic mean free paths based on 
%       calculations for 31 materials. Surface Science Letters, 192(1), 
%       L849–L857. https://doi.org/10.1016/0167-2584(87)90829-2
%
%   IN:
%   -   ke_dat:     N×1 vector of the input electron kinetic energy (for PES; KE = HV - BE - PHI) [eV]
%   -   material:  	string of the material whose imfp is to be determined; e.g. "Si", "SiO2", "Al2O3"...
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]

%% Default parameters (Parameters for Silicon)
if nargin < 2; material = "Si";  end
if isempty(material); material = "Si"; end
%% - 1 - Extracting the material parameters from the materials database
% - Ensure the input is a string
material = string(material);
% - Extracting the material properties
material_props = get_mpd_props(material);
% - Extracting the material properties required to determine the TPP-2M IMFP
rho     = material_props.DENSITY;
Nv      = material_props.ELECT_VALENCY;
M       = material_props.ATOM_MASS;
Egap    = material_props.ELE_BGAP;
%% - 2 - Determination of the IMFP TPP-2M
imfp = eimfp_tpp2m(ke_dat, rho, Nv, M, Egap);
%% Ensuring the imfp is a column vector
if size(imfp, 2) >1; imfp = imfp'; end
end