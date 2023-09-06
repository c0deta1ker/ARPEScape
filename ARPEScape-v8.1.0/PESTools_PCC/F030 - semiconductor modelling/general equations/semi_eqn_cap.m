function capacitance = semi_eqn_cap(area, epsR, distance)
% capacitance = semi_eqn_cap(area, epsR, distance)
%   Function that determines the capacitance of two parallel plates, with a
%   dielectric material placed between them. The thickness and relative
%   permittivity of the dieelctric can be defined by the user, as well as
%   the cross sectional area of the capacitor plates.
%
%   REQ. FUNCTIONS: 
%   - physics_constants();
%
%   IN:
%   -   capacitance:    capacitance [Farads]
%
%   OUT:
%   -   area:           cross-sectional area between the capacitor plates [m^2]
%   -   epsR:           relative permittivity of dielectric.
%   -   distance:       distance between the plates of the capacitor [m].

%% Defining physical constants
pc = physics_constants();
%% - 1 - Defining the capacitance equation
capacitance = pc.eps0*epsR*area/distance;
end