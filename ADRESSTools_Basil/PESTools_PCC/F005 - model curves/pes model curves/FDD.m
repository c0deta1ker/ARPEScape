function int = FDD(xdat, fdd_ef, fdd_T)
% int = FDD(xdat, fdd_ef, fdd_T)
%   This function determines the Fermi-Dirac Distribution (FDD) given the 
%   fermi-energy (fdd_ef) and temperature (fdd_T). This is the pure, theoretical
%   FDD function.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       Nx1 column vector of the input x-axis domain (binding energy for PES)
%   -   fdd_ef:   	scalar of the Fermi-level position (eV).
%   -   fdd_T:     	scalar of the temperature (K).
%
%   OUT:
%   -   int:        Nx1 column vector of the output FDD curve profile.

%% Default parameters
if nargin < 3; fdd_T = 12;  end
if nargin < 2; fdd_ef = 0; fdd_T = 12; end
if isempty(fdd_ef); fdd_ef = 0; end
if isempty(fdd_T); fdd_T = 12; end
% - Validity check on the input parameters
if fdd_T < 0; fdd_T = 0; end

%% 1 - Determination of the Fermi-Dirac Distribution
% Ensuring xdat is a column vector
if size(xdat, 2) >1; xdat = xdat'; end
% Evaluating the FDD
int = 1 ./ (exp((xdat - fdd_ef) ./ (8.617e-5 .* fdd_T)) + 1);
% If isnan, return zero
int(isnan(int)) = 0;
end