function imfp = eimfp_tpp2m(ke_dat, rho, Nv, M, Egap)
% imfp = eimfp_tpp2m(ke_dat, rho, Nv, M, Egap)
%   Function that determines the electron inelastic mean free path (IMFP) 
%   based on the TPP-2M energy compensation factor [1], which is a more
%   accurate method than the standard universal curve determination. The
%   input parameters of the materials are now required. The TPP-2M equation
%   is also relativistically corrected here, as per reference [2].
%
%   REFERENCE:
%   [1] Tanuma, S., Powell, C. J., & Penn, D. R. (1987). 
%       Proposed formula for electron inelastic mean free paths based on 
%       calculations for 31 materials. Surface Science Letters, 192(1), 
%       L849–L857. https://doi.org/10.1016/0167-2584(87)90829-2
%   [2] Shinotsuka H, Tanuma S, Powell CJ, Penn DR. Calculations of 
%       electron inelastic mean free paths. XII. Data for 42 inorganic 
%       compounds over the 50 eV to 200 keV range with the full Penn 
%       algorithm. Surf Interface Anal. 2018;51(4):427-457. 
%       doi:10.1002/sia.6598
%   [1] M. P. Seah, “An accurate and simple universal curve for the 
%       energy-dependent electron inelastic mean free path,” 
%       Surf. Interface Anal., vol. 44, no. 4, pp. 497–503, 2012, 
%       doi: 10.1002/sia.4816.

%
%   IN:
%   -   ke_dat:     N×1 vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%   -   rho:        scalar of the density of the material (g/cc)
%   -   Nv:         scalar of the number of valence electrons per atom (for an element)
%   -   M:          scalar of the atomic or molecular weight (in amu == g/mol)
%   -   Egap:       scalar of the band gap energy (eV)
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]

%% Default parameters (Parameters for Silicon)
if nargin < 5; Egap = 1.12000; end
if nargin < 4; M = 28.085000;   end
if nargin < 3; Nv = 4;  end
if nargin < 2; rho = 2.3300000; end
if isempty(Egap);   Egap = 1.12000; end     % eV
if isempty(M);      M = 28.085000; end      % amu == g/mol
if isempty(Nv);     Nv = 4; end             % integer
if isempty(rho); 	rho = 2.3300000; end    % g/cc
%% - 1 - Determination of the IMFP TPP-2M
% If the kinetic energy is negative, assume it is zero
ke_dat(ke_dat<0) = 0;
% Evaluating the imfp
Ep          = 28.816 .* sqrt(Nv .* rho ./ M); % eV
alpha       = (1 + ke_dat ./ 1021999.8) ./ (1 + ke_dat ./ 510998.9).^2;     % relativistic correction
beta        = -1.0 + 9.44 ./ sqrt(Ep.^2 + Egap.^2) + 0.69 .* rho.^(0.1);    % (eV−1 nm−1)
gamma       = 0.191 .* rho.^(-0.5);     % eV−1
U       	= (Ep ./ 28.816).^2; 
D           = 534 - 208 .* U;          % eV nm−1
C           = 19.7 - 9.1 .* U;         % nm−1
% Relativistic corrected equation
imfp = alpha .* ke_dat ./ (Ep.^2 .* (beta .* log(alpha .* gamma .* ke_dat) - (C./ke_dat) + (D./(ke_dat.^2)))); % IMFP in nanometres
% If isnan, return zero
imfp(isnan(imfp)) = 0;
% Convert IMFP from nm to Angstroms
imfp = 10 .* imfp;
%% Ensuring the imfp is a column vector
if size(imfp, 2) >1; imfp = imfp'; end
end