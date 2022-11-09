function kf = moscap_eqn_kf(ne, ni)
% kf = moscap_eqn_kf(ne, ni)
%
%   REQ. FUNCTIONS: 
%   - physics_constants();
%
%   IN:
%   -   kf:             Fermi wave-vector of the conduction electrons [Ang^-1]
%
%   OUT:
%   -   ne:             charge accumulated on the capacitor plates [cm^-2]
%   -   ni:             initial charge on the capacitor plates [cm^-2]
%% Default parameters
if nargin < 2; ni = 0; end
if isempty(ni); ni = 0; end
%% Defining physical constants
pc = physics_constants();
%% - 1 - Defining the Fermi wave-vector
kf = sqrt(2*pi*(ne + ni)*10000)*1e-10;
end