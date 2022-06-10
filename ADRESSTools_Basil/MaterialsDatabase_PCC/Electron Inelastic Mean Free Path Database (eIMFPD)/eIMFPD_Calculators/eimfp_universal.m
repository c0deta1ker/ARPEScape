function imfp = eimfp_universal(ke_dat)
% imfp = eimfp_universal(ke_dat)
%   Function that determines the universal electron inelastic mean free
%   path (IMFP) in elements based on the KE^0.5 equation. This gives a
%   good, first order approximation of the IMFP values, with the only
%   required input being the electron kinetic energy.
%
%   IN:
%   -   ke_dat:  	N×1 column vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]

%% - 1 - Determination of the IMFP using universal formula
% If the kinetic energy is negative, assume it is zero
ke_dat(ke_dat<0) = 0;
% Evaluating the imfp
imfp = (143 ./ ke_dat.^2) + 0.054 .* sqrt(ke_dat);       % IMFP in nano-metres
% If isnan, return zero
imfp(isnan(imfp)) = 0;
% Convert IMFP from nm to Angstroms
imfp = 10 .* imfp;
%% Ensuring the imfp is a column vector
if size(imfp, 2) >1; imfp = imfp'; end
end