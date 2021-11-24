function k = k_calc(eb, k0, eb0, mstar)
% k = k_calc(eb, k0, eb0, mstar)
%   Function that evaluates the wave-vector of a parabolic dispersion;
%   k0, eb0 are the initial values at the base of the parabola and mstar
%   defines the curvature of the parabolic disperison.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   eb:        Nx1 vector that yields the binding energy values [eV].
%   -   k0:        scalar of the initial wave-vector at eb = 0 eV [Ang^-1].
%   -   eb0:       scalar of the initial binding energy at k = 0 Ang^-1 [eV].
%   -   mstar:     scalar of the effective mass which defines the curvature of the parabola [mstar*me].
%
%   OUT:
%   -   k:         Nx1 vector of the output wave-vector, k, values [Ang^-1].

%% Default parameters
if nargin < 4; mstar = 1; end
if nargin < 3; eb0 = 0; end
if nargin < 2; k0 = 0; end
if isempty(mstar); mstar = 1; end
if isempty(eb0); eb0 = 0; end
if isempty(k0); k0 = 0; end

%% - 1 - Determination of the output wave-vectors
% Ensuring eb is a column vector
if size(eb, 2) > 1 && size(eb, 1) == 1; eb = eb'; end
% Defining the physical constants
hbar    = 4.135e-15;    % Units of [eV s^−1]
% Calculating the wave-vector values
k       = sqrt((mstar / 3.8180) .* (eb - eb0)) + k0;
% If isnan, return zero
k(isnan(k)) = 0;
end