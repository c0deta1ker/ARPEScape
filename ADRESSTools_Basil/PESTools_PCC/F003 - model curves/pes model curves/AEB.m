function int = AEB(xdat, x0, asym)
% int = AEB(xdat, x0, asym)
%   Function that evaluates an asymmetric exponential blend (AEB), that is 
%   used to create the asymmetric Voigt-type lineshapes. The asymmetric profile
%   obtained from the blend is: Y(x) = GL(x) + AEB(x) * (1 - GL(x)).
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   xdat:    	Nx1 column vector of the input domain (binding energy for XPS)
%   -   x0:         scalar of the position where the AEB will start.
%   -   asym:       asymmetry ratio; 0 for none, 1 is for maximum
%
%   OUT:
%   -   int:        Nx1 column vector of the asymmetric exponential blend intensity

%% Default parameters
if nargin < 3; asym = 0.5; end
if nargin < 2; x0 = 0; asym = 0.5; end
if isempty(x0); x0 = 0; end
if isempty(asym); asym = 0.5; end

%% - 1 - Determination of the Exponential Asymmetric Blend intensities
% Ensuring xdat is a column vector
if size(xdat, 2) >1; xdat = xdat'; end
% Evaluating the asymmetric blend
K       = 1;    % Coefficient
int     = exp(-K/asym.*abs(xdat - x0));
% For all values where xdat > x0, set to zero
[~, indx]       = min(abs(xdat - x0));
int(indx:end)   = 0;
end