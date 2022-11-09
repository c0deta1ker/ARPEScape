function current = semi_eqn_leak(voltage, area, distance, resistivity)
% current = semi_eqn_leak(voltage, area, distance, resistivity)
%   Function that determines the leakage current between two parallel
%   plates of a capacitor, through the intermediate dieelctric layer.
%
%   REQ. FUNCTIONS: 
%   - physics_constants();
%
%   IN:
%   -   current:      	leakage current flowing through capacitor [A]
%
%   OUT:
%   -   voltage:        potential difference across the capacitor [V].
%   -   area:           cross-sectional area of the capacitor [m^2].
%   -   distance:       distance between the plates of the capacitor [m].
%   -   resistivity:    resistivity of the dielectric in the capacitor [Ohm m]
%% Defining physical constants
pc = physics_constants();
%% - 1 - Defining the leakage current equation
current = voltage.*area./resistivity./distance;
end