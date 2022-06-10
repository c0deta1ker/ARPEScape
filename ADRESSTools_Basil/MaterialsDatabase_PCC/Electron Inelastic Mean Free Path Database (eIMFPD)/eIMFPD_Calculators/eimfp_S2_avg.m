function imfp = eimfp_S2_avg(ke_dat, Z)
% imfp = eimfp_S2_avg(ke_dat, Z)
%   Function that determines the electron inelastic mean free path (IMFP) 
%   based on the S2 equation described by M. P. Seah [1]. The IMFP here
%   only depends on the value of Z.
%
%   REFERENCE:
%   [1] M. P. Seah, “An accurate and simple universal curve for the 
%       energy-dependent electron inelastic mean free path,” 
%       Surf. Interface Anal., vol. 44, no. 4, pp. 497–503, 2012, 
%       doi: 10.1002/sia.4816.
%
%   IN:
%   -   ke_dat:     N×1 vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%   -   Z:          scalar of the atomic mass number (Z) (or average for compounds)
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]

%% Default parameters (Parameters for Silicon)
if nargin < 2;  Z = 14; end
if isempty(Z); 	Z = 14; end

%% - 1 - Determination of the IMFP TPP-2M
% If the kinetic energy is negative, assume it is zero
ke_dat(ke_dat<0) = 0;
% Evaluating the imfp
a           = 2.5;  % Average atomic spacing in Angstroms
% Relativistic corrected equation
imfp = a .* ((1.52 + 0.167 .* Z.^(0.5) + 0.0394 .* ke_dat.^(0.872)) ./ Z.^(0.3)); % IMFP in Angstroms
% If isnan, return zero
imfp(isnan(imfp)) = 0;
%% Ensuring the imfp is a column vector
if size(imfp, 2) >1; imfp = imfp'; end
end