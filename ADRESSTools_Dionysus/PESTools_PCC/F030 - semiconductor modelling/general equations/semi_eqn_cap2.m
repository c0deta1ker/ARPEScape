function capacitance = semi_eqn_cap2(area, epsR_1, distance_1, epsR_2, distance_2)
% capacitance = semi_eqn_cap2(area, epsR_1, distance_1, epsR_2, distance_2)
%   Function that determines the capacitance for two dielectric materials
%   that are placed in series, between the capacitor plates. The user can
%   define the relative perimittivity and thickness of each layer. The
%   cross sectional area of the capacitor is constant for the two layers.
%
%   REQ. FUNCTIONS: 
%   - physics_constants();
%
%   IN:
%   -   capacitance:    capacitance [Farads]
%
%   OUT:
%   -   area:           cross-sectional area between the capacitor plates [m^2]
%   -   epsR_1:         relative permittivity of dielectric layer 1.
%   -   distance_1:     thickness of layer 1 between the capacitor plates [m].
%   -   epsR_2:         relative permittivity of dielectric layer 2.
%   -   distance_2:     thickness of layer 2 between the capacitor plates [m].

%% Defining physical constants
pc = physics_constants();
%% - 1 - Defining the capacitance equation
capacitance = ((epsR_1*epsR_2)/(epsR_1*distance_2+epsR_2*distance_1))*(pc.eps0*area);
end