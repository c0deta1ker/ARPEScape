function ydat = FDDG(xdat, fdd_ef, fdd_T, fdd_fwhm)
% ydat = FDDG(xdat, fdd_ef, fdd_T, fdd_fwhm)
%   This function determines the Fermi-Dirac Distribution (FDD) given the 
%   fermi-energy (fdd_ef) and temperature (fdd_T). A convolution is then made with
%   a Gaussian, with a well-defined full-width half-maximum (fdd_fwhm) to
%   include the broadening introducted by the beamline+analyzer resolution.
%
%   REQ. FUNCTIONS:
%   -   FDD()
%   -   Gauss()
%
%   IN:
%   -   xdat:       N×1 column vector of the input domain (binding energy for PES).
%   -   fdd_ef:   	scalar of the Fermi-level position (eV).
%   -   fdd_T:     	scalar of the temperature (K).
%   -   fdd_fwhm: 	scalar of the full-width at half-maximum (FWHM) of the Gaussian broadening term (eV).
%
%   OUT:
%   -   ydat:       N×1 column vector of the intensity range (intensity for PES).

%% Default parameters
if nargin < 2; fdd_ef = 0; end
if nargin < 3; fdd_T = 12;  end
if nargin < 4; fdd_fwhm = 0;  end
if isempty(fdd_ef); fdd_ef = 0; end
if isempty(fdd_T); fdd_T = 12; end
if isempty(fdd_fwhm); fdd_fwhm = 0; end
%% Validity checks on the input parameters
if fdd_T < 0; fdd_T = 0; end
if fdd_fwhm < 0; fdd_fwhm = 0; end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end

%% 1 - Determining the FDD and Gaussian functions to use
% - Defining a large domain that spans the desired domain
xdat_large = linspace(min(xdat(:)) - 2, max(xdat(:)) + 2, 1e4)';
% - FDD function
int_fdd     = FDD(xdat_large, fdd_ef, fdd_T);
% - Gaussian function
int_Gauss  	= Gauss(xdat_large, fdd_ef, 1, fdd_fwhm);
%% 2 - Convolving the FDD and Gaussian functions
y_conv      = conv(int_fdd, int_Gauss, 'same'); 
y_conv      = y_conv / max(y_conv);
% - Eliminate edge effects of the convolution
for i = 0:size(y_conv, 1)-1
    y_conv_val = y_conv(end-i); if y_conv_val == 1; y_conv(1:end-i) = 1; break; end
end
% - Validity check that the Fermi-level is the same
[~, indx]       = min(abs(y_conv - 0.5));
ef_shift        = xdat_large(indx) - fdd_ef;
% - Assigning the best fit Fermi function
fdd_fit         = fit(xdat_large-ef_shift, y_conv, 'linearinterp');  
%% 3 - Defining the final intensity array
ydat = fdd_fit(xdat);
% If isnan, return zero
ydat(isnan(ydat)) = 0;

end