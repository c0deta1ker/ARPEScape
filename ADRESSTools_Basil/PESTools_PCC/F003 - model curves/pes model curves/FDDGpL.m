function ydat = FDDGpL(xdat, fdd_ef, fdd_T, fdd_fwhm, lin_grad, lin_inter, const_bgrnd)
% ydat = FDDGpL(xdat, fdd_ef, fdd_T, fdd_fwhm, lin_grad, lin_inter, const_bgrnd)
%   This function determines the Fermi-Dirac Distribution (FDD) given the 
%   fermi-energy (fdd_ef) and temperature (fdd_T). A convolution is then made with
%   a Gaussian, with a well-defined full-width half-maximum (fdd_fwhm) to
%   include the broadening introducted by the beamline+analyzer resolution.
%   Furthermore, a linear background is MULTIPLIED to accound for any
%   linear distortion of the Fermi-edge, defined by the gradient and intercept.
%   An additional constant background is also added (const_bgrnd).
%
%   REQ. FUNCTIONS: 
%   -   FDDG()
%
%   IN:
%   -   xdat:           N×1 column vector of the input domain (binding energy for PES).
%   -   fdd_ef:         scalar of the Fermi-level position (eV).
%   -   fdd_T:          scalar of the temperature (K).
%   -   fdd_fwhm:       scalar of the full-width at half-maximum (FWHM) of the Gaussian broadening term (eV).
%   -   lin_grad:       scalar of the gradient of the linear background.
%   -   lin_offset:  	scalar of the y-intercept of the linear background.
%   -   const_bgrnd:  	scalar of the constant background y-offset value.
%
%   OUT:
%   -   ydat:           N×1 column vector of the intensity range (intensity for PES).

%% Default parameters
if nargin < 2; fdd_ef = 0; end
if nargin < 3; fdd_T = 12; end
if nargin < 4; fdd_fwhm = 0; end
if nargin < 5; lin_grad = 0; end
if nargin < 6; lin_inter = 1; end
if nargin < 7; const_bgrnd = 0; end
if isempty(fdd_ef);     fdd_ef = 0; end
if isempty(fdd_T);      fdd_T = 12; end
if isempty(fdd_fwhm);   fdd_fwhm = 0; end
if isempty(lin_grad);   lin_grad = 0; end
if isempty(lin_inter);  lin_inter = 1; end
if isempty(const_bgrnd);const_bgrnd = 0; end
%% Validity checks on the input parameters
if fdd_T < 0; fdd_T = 0; end
if fdd_fwhm < 0; fdd_fwhm = 0; end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end

%% 1 - Determining the FDD, linear and background intensities
% - FDD function
int_fdd     = FDDG(xdat, fdd_ef, fdd_T, fdd_fwhm);
% - Linear function
int_linear	= lin_grad .* xdat + lin_inter;
% - Background value
int_bgrnd  	= const_bgrnd;
%% 2 - Combining all the functions
ydat        = int_bgrnd + int_linear .* int_fdd;
% If isnan, return zero
ydat(isnan(ydat)) = 0;

end