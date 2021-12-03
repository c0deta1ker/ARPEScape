function eb = eb_calc(k, k0, eb0, mstar)
% eb = eb_calc(k, k0, eb0, mstar)
%   Function that evaluates the binding energy of a parabolic dispersion;
%   k0, eb0 are the initial values at the base of the parabola and mstar
%   defines the curvature of the parabolic disperison.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   k:         Nx1 vector that yields the wave-vector values [Ang^-1].
%   -   k0:        scalar of the initial wave-vector at eb = 0 eV [Ang^-1].
%   -   eb0:       scalar of the initial binding energy at k = 0 Ang^-1 [eV].
%   -   mstar:     scalar of the effective mass which defines the curvature of the parabola [mstar*me].
%
%   OUT:
%   -   eb:        Nx1 vector that yields the binding energy values [eV].

%% Default parameters
if nargin < 4; mstar = 1; end
if nargin < 3; eb0 = 0; end
if nargin < 2; k0 = 0; end
if isempty(mstar); mstar = 1; end
if isempty(eb0); eb0 = 0; end
if isempty(k0); k0 = 0; end

%% - 1 - Determination of the output wave-vectors
% Ensuring eb is a column vector
if size(k, 2) > 1 && size(k, 1) == 1; k = k'; end
% Defining the physical constants
hbar    = 4.135e-15;    % Units of [eV s^−1]
% Calculating the binding energy values
eb = (3.8180 / mstar) * (k - k0).^2 + eb0;
% If isnan, return zero
eb(isnan(eb)) = 0;
end