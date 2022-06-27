function ydat = FDD(xdat, fdd_ef, fdd_T)
% ydat = FDD(xdat, fdd_ef, fdd_T)
%   This function determines the Fermi-Dirac Distribution (FDD) given the 
%   fermi-energy (fdd_ef) and temperature (fdd_T). This is the pure, theoretical
%   FDD function.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       N×1 column vector of the input domain (binding energy for PES).
%   -   fdd_ef:   	scalar of the Fermi-level position (eV).
%   -   fdd_T:     	scalar of the temperature (K).
%
%   OUT:
%   -   ydat:       N×1 column vector of the intensity range (intensity for PES).

%% Default parameters
if nargin < 2; fdd_ef = 0; end
if nargin < 3; fdd_T = 12;  end
if isempty(fdd_ef); fdd_ef = 0; end
if isempty(fdd_T); fdd_T = 12; end
%% Validity checks on the input parameters
if fdd_T < 0; fdd_T = 0; end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end

%% 1 - Determination of the Fermi-Dirac Distribution
% Evaluating the FDD
ydat = 1 ./ (exp((xdat - fdd_ef) ./ (8.617e-5 .* fdd_T)) + 1);
% If isnan, return zero
ydat(isnan(ydat)) = 0;

end