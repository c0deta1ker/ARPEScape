function imfp = eimfp_S1(ke_dat, rho, M, Egap, Z)
% imfp = eimfp_S1(ke_dat, rho, M, Egap, Z)
%   Function that determines the electron inelastic mean free path (IMFP) 
%   based on the S1 equation described by M. P. Seah [1]. The IMFP here
%   additionally depends on the value of Z and this formalism is compatible
%   for elemental or binary materials.
%
%   REFERENCE:
%   [1] M. P. Seah, “An accurate and simple universal curve for the 
%       energy-dependent electron inelastic mean free path,” 
%       Surf. Interface Anal., vol. 44, no. 4, pp. 497–503, 2012, 
%       doi: 10.1002/sia.4816.
%   
%   IN:
%   -   ke_dat:     N×1 vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%   -   rho:        scalar of the density of the material (g/cc)
%   -   M:          scalar of the atomic or molecular weight (in amu == g/mol)
%   -   Egap:       scalar of the band gap energy (eV)
%   -   Z:          scalar of the atomic mass number of the element (or average for compound) (Z) 
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]

%% Default parameters (Parameters for Silicon)
if nargin < 5; Z = 14; end
if nargin < 4; Egap = 1.12000; end
if nargin < 3; M = 28.085000;   end
if nargin < 2; rho = 2.3300000; end
if isempty(Z);      Z = 14; end
if isempty(Egap);   Egap = 1.12000; end     % eV
if isempty(M);      M = 28.085000; end      % amu == g/mol
if isempty(rho); 	rho = 2.3300000; end    % g/cc
% -- Defining the constants
NA  = 6.02214086e23;    % Avogadro constant
%% - 1 - Determination of the IMFP using S1 formalism (M P Seah)
% If the kinetic energy is negative, assume it is zero
ke_dat(ke_dat<0) = 0;
% Evaluating the imfp
a           = ((1e21 .* M) ./ (rho .* NA)).^(1./3); % nm
% Relativistic corrected equation
imfp = (4 + 0.44 .* Z.^(0.5) + 0.104 .* ke_dat.^(0.872)) .* a.^(1.7) ./ (Z.^(0.3) .* (1 - 0.02 .*Egap)); % IMFP in nanometres
% If isnan, return zero
imfp(isnan(imfp)) = 0;
% Convert IMFP from nm to Angstroms
imfp = 10 .* imfp;
%% Ensuring the imfp is a column vector
if size(imfp, 2) >1; imfp = imfp'; end
end