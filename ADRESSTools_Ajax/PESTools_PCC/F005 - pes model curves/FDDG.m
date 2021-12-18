function int = FDDG(xdat, fdd_ef, fdd_T, fdd_fwhm)
% int = FDDG(xdat, fdd_ef, fdd_T, fdd_fwhm)
%   This function determines the Fermi-Dirac Distribution (FDD) given the 
%   fermi-energy (fdd_ef) and temperature (fdd_T). A convolution is then made with
%   a Gaussian, with a well-defined full-width half-maximum (fdd_fwhm) to
%   include the broadening introducted by the beamline+analyzer resolution.
%
%   REQ. FUNCTIONS:
%   -   FDD(xdat, ef, T)
%   -   Gauss(xdat, x0, peak, fwhm)
%
%   IN:
%   -   xdat:       Nx1 column vector of the input domain (binding energy for XPS)
%   -   fdd_ef:   	scalar of the Fermi-level position (eV).
%   -   fdd_T:    	scalar of the temperature (K).
%   -   fdd_fwhm: 	scalar of the full-width at half-maximum (FWHM) of the Gaussian broadening term (eV).
%
%   OUT:
%   -   int:        Nx1 column vector of the output FDD curve profile.

%% Default parameters
if nargin < 4; fdd_fwhm = 0;  end
if nargin < 3; fdd_T = 12; fdd_fwhm = 0;  end
if nargin < 2; fdd_ef = 0; fdd_T = 12; fdd_fwhm = 0; end
if isempty(fdd_fwhm); fdd_fwhm = 0; end
if isempty(fdd_T); fdd_T = 12; end
if isempty(fdd_ef); fdd_ef = 0; end
% - Validity check on the input parameters
if fdd_T < 0; fdd_T = 0; end
if fdd_fwhm < 0; fdd_fwhm = 0; end

%% 1 - Determining the FDD and Gaussian functions to use
% Ensuring xdat is a column vector
if size(xdat, 2) >1; xdat = xdat'; end
% - Defining a large domain that spains the desired domain
xdat_large = linspace(min(xdat(:)) - 5, max(xdat(:)) + 5, 1e3)';
% - FDD function
int_fdd     = FDD(xdat_large, fdd_ef, fdd_T);
% - Gaussian function
int_Gauss  	= Gauss(xdat_large, fdd_ef, 1, fdd_fwhm);
%% 2 - Convolving the FDD and Gaussian functions
y_conv      = conv(int_Gauss, int_fdd, 'same'); 
y_conv      = y_conv / max(y_conv);
% - Eliminate edge effects of the convolution
for i = 1:size(y_conv, 2)
    y_conv_val = y_conv(i); if y_conv_val == 1; y_conv(1,1:i) = 1; break; end
end
% - Assigning the best fit Fermi function
fdd_fit = fit(xdat_large, y_conv, 'linearinterp');  
%% 3 - Defining the final intensity array
int = fdd_fit(xdat);
% If isnan, return zero
int(isnan(int)) = 0;
end