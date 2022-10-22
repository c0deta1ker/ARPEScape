function ydat = AEB(xdat, x0, asym)
% ydat = AEB(xdat, x0, asym)
%   Function that evaluates an asymmetric exponential blend (AEB), that is 
%   used to create the asymmetric Voigt-type lineshapes. The asymmetric profile
%   obtained from the blend is: Y(x) = GL(x) + AEB(x) * (1 - GL(x)).
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:       N×1 column vector of the input domain (binding energy for PES).
%   -   x0:         scalar of the position where the AEB will start.
%   -   asym:       asymmetry ratio; 0 for none, 0->1 is for LHS asymmetry, -1->0 is for RHS asymmetry.
%
%   OUT:
%   -   ydat:       N×1 column vector of the intensity range (intensity for PES).

%% Default parameters
if nargin < 2;      x0      = 0.0; end
if nargin < 3;      asym    = 0.5; end
if isempty(x0);     x0      = 0; end
if isempty(asym);   asym    = 0.5; end
%% Validity checks on the input parameters
% -- Padding the asymmetry to be between -1 -> 1
if asym < -1; asym = -1; end
if asym > 1; asym = 1; end
%% Ensuring the data is in the form of 1D column vectors
if size(xdat, 2) > 1; xdat = xdat'; end

%% - 1 - Determination of the Exponential Asymmetric Blend intensities
% Evaluating the asymmetric blend
if asym > 0 || asym == 0
    K       = 1;    % Coefficient
    ydat    = exp(-K/asym.*abs(xdat - x0));
    % For all values where xdat > x0, set to zero
    [~, indx]       = min(abs(xdat - x0));
    ydat(indx:end)	= 0;
elseif asym < 0
    K       = 1;    % Coefficient
    ydat    = exp(K/asym.*abs(xdat - x0));
    % For all values where xdat < x0, set to zero
    [~, indx]       = min(abs(xdat - x0));
    ydat(1:indx)	= 0;
end
% If isnan, return zero
ydat(isnan(ydat)) = 0;

end