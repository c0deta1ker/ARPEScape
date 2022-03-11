function zdat = ARPESCurveNorm(xdat, ydat, cTYPE, INT, XLOC, YLOC, XFWHM, YFWHM, MSTAR)
% zdat = ARPESCurve(xdat, ydat, cTYPE, INT, XLOC, YLOC, XFWHM, YFWHM, MSTAR)
%   Function that evaluates a generic angle-resolved photoelectron 
%   spectroscopy (ARPES) curve in 2D. This can be used to model an ARPES
%   spectrum, that is built up by a sum of individual Gaussians.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:    	1xM row vector of the x-axis domain
%   -   ydat:    	Nx1 column vector of the y-axis domain
%   -   cTYPE:      type of curve to use for fitting. Default: "G2DA" ("G2D", "L2D")
%   -   INT:    	scalar of the peak intensity of 2D PE curve.
%   -   XLOC:      	scalar of the x-location of the min/max of the 2D parabolic ARPES dispersion [Ang^-1].
%   -   YLOC:      	scalar of the y-location of the min/max of the 2D parabolic ARPES dispersion [eV].
%   -   XFWHM:     	scalar of the x-axis FWHM for each Gaussian (k-resolution) [Ang^-1]
%   -   YFWHM:     	scalar of the y-axis FWHM for each Gaussian (Eb-resolution) [eV]
%   -   MSTAR:     	scalar of the effective mass, which governs the curvature of the parabola.
%
%   OUT:
%   -   zdat:     	N x M matrix of the output 2D parabolic ARPES dispersion.

%% Default parameters
% Default based on inputs
if nargin < 9; MSTAR = 0.20;  end
if nargin < 8; YFWHM = 0.25; end
if nargin < 7; XFWHM = 0.02; end
if nargin < 6; YLOC	= 0.00; end
if nargin < 5; XLOC	= 0.00; end
if nargin < 4; INT  = 1.00; end
if nargin < 3; cTYPE = "G2DA"; end
% Default based on empty inputs
if isempty(cTYPE);  cTYPE   = "G2DA"; end
if isempty(INT);    INT     = 1.00; end
if isempty(XLOC);   XLOC    = 0.00; end
if isempty(YLOC);   YLOC    = 0.00; end
if isempty(XFWHM);  XFWHM   = 0.02; end
if isempty(YFWHM);  YFWHM   = 0.25; end
if isempty(MSTAR);  MSTAR   = 0.20; end

%% - 1 - Determination of the angle-resolved photoemission spectrum
% Ensuring xdat is a row vector
if size(xdat, 1) > 1; xdat = xdat'; end
% Ensuring ydat is a column vector
if size(ydat, 2) > 1; ydat = ydat'; end
% Filing through all the curves
for j = 1:length(cTYPE)
    % Initialising the data
    k_val   = xdat;
    eb_val  = eb_calc(xdat, XLOC(j), YLOC(j), MSTAR(j));
    P_curve = zeros(length(ydat), length(xdat), length(k_val));
    % Extracting the primary curve components
    for i = 1:length(k_val)
        if cTYPE == "G2D"
            P_curve(:,:,i)  = G2D(xdat, ydat, INT(j), k_val(i), eb_val(i), XFWHM(j)); 
        elseif cTYPE == "L2D"
            P_curve(:,:,i)  = L2D(xdat, ydat, INT(j), k_val(i), eb_val(i), XFWHM(j)); 
        elseif cTYPE == "G2DA"
            P_curve(:,:,i)  = G2DA(xdat, ydat, INT(j), k_val(i), eb_val(i), XFWHM(j), YFWHM(j)); 
        end
    end
    % Summing the components together
    zdat_temp{j} = sum(P_curve, 3);
    % If isnan, return zero
    zdat_temp{j}(isnan(zdat_temp{j})) = 0;
end
% Summing over all the data
zdat = zdat_temp{1};
if length(zdat_temp) > 1; for i = 2:length(zdat_temp); zdat = zdat + zdat_temp{i}; end; end

%% - 2 - Normalising across all of the EDCs
for i = 1:size(zdat,2)
    zdat(:,i) = zdat(:,i) ./ max(zdat(:,i));
end

end