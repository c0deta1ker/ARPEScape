function ne = moscap_eqn_ne(voltage, epsR, distance)
% ne = moscap_eqn_ne(voltage, epsR, distance)
%
%   REQ. FUNCTIONS: 
%   - physics_constants();
%
%   IN:
%   -   ne:             charge accumulated on the capacitor plates [m^-2]
%
%   OUT:
%   -   voltage:        potential difference across the capacitor [V].
%   -   epsR:           relative permittivity of dielectric.
%   -   distance:       distance between the plates of the capacitor [m].
%% Defining physical constants
pc = physics_constants();
%% - 1 - Defining the charge accumulation equation
ne = pc.eps0*epsR*voltage/distance/pc.e;
end