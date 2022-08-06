function capacitance = moscap_eqn_cap(epsR, area, distance)
% capacitance = moscap_eqn_cap(voltage, eps_R, distance)
%
%   REQ. FUNCTIONS:
%
%   IN:
%   -   capacitance:    capacitance [Farads]
%
%   OUT:
%   -   voltage:        potential difference across the capacitor [V].
%   -   epsR:           relative permittivity of dielectric.
%   -   distance:       distance between the plates of the capacitor [m].
%% Defining physical constants
pc = physics_constants();
%% - 1 - Defining the capacitance equation
capacitance = pc.eps0*epsR*area/distance;
end