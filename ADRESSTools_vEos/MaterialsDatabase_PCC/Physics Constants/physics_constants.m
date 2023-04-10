function pc = physics_constants()
% pc = phys_consts()
%   This function creates a MATLAB structure that returns all the physical
%   constants necessary in Physics calculations. The physical constants are
%   defined using the International System of Units (SI).
%
%   IN: (none)
%
%   OUT:
%   -   pc:     MATLAB structure that contains all of the relevant physics constants.
%% - 1 - DEFINING ALL THE PHYSICS CONSTANTS
pc          = struct();
pc.G        = 6.6743015e-11;        % Newtonian constant of gravitation [m3⋅kg−1⋅s−2]
pc.c        = 299792458.;           % speed of light in vacuum [m⋅s−1]
pc.h        = 6.62607015e-34;       % Planck constant [J⋅Hz−1]
pc.hbar     = pc.h/2/pi;            % reduced Planck constant [J⋅s]
pc.mu0      = 1.25663706212e-6;     % vacuum magnetic permeability [N⋅A−2]
pc.Z0       = pc.mu0*pc.c;          % characteristic impedance of vacuum [Ω]
pc.eps0     = 1/pc.mu0/pc.c^2;      % vacuum electric permittivity [F⋅m−1]
pc.ke       = 1/4/pi/pc.eps0;       % Coulomb constant [N⋅m2⋅C−2]
pc.kB       = 1.380649e-23;         % Boltzmann constant [J⋅K−1]
pc.sigma    = 5.670374419e-8;       % Stefan–Boltzmann constant [W⋅m−2⋅K−4]
pc.c1       = 3.741771852e-16;      % first radiation constant [W⋅m2]
pc.c1L      = 1.191042972e-16;      % first radiation constant for spectral radiance [W⋅m2⋅sr−1]
pc.c2       = pc.h*pc.c/pc.kB;      % second radiation constant [m⋅K]
pc.b        = 2.897771955e-3;       % Wien wavelength displacement law constant [m⋅K]
pc.bprime   = pc.c/pc.b;            % Wien frequency displacement law constant [Hz⋅K−1]
pc.bentropy = 3.002916077e-3;       % Wien entropy displacement law constant [m⋅K]
pc.e        = 1.602176634e-19;      % elementary charge [C]
pc.G0       = 2*pc.e^2/pc.h;        % conductance quantum [S]
pc.G0i      = pc.h/2/pc.e^2;        % inverse conductance quantum [Ω]
pc.Rk       = pc.h/pc.e^2;          % von Klitzing constant [Ω]
pc.Kj       = 2*pc.e/pc.h;          % Josephson constant [Hz⋅V−1]
pc.Phi0     = pc.h/2/pc.e;          % magnetic flux quantum [Wb]
pc.alpha    = pc.e^2/4/pi/pc.eps0/pc.hbar/pc.c;     % fine-structure constant
pc.me       = 9.109383701528e-31;   % electron mass [kg]
pc.mp       = 1.6726219236951e-27;  % proton mass [kg]
pc.mn       = 1.6749274980495e-27;  % neutron mass [kg]
pc.mmu      = 1.88353162742e-28;    % muon mass [kg]
pc.mtau     = 3.1675421e-27;        % tau mass [kg]
pc.mt       = 3.078453e-25;         % top quark mass [kg]
pc.mpe      = 1836.15267343;        % proton-to-electron mass ratio
pc.mWZ      = 0.8815317;            % W-to-Z mass ratio
pc.ge       = -2.0023193043625635;  % electron g-factor
pc.gmu      = -2.002331841813;      % muon g-factor
pc.gp       = 5.585694689316;       % proton g-factor
pc.qc       = pc.h/2/pc.me;         % quantum of circulation [m2⋅s−1]
pc.muB      = 9.2740100783e-24;     % Bohr magneton
pc.muN      = 5.0507837461e-27;     % nuclear magneton
pc.re       = pc.e^2*pc.ke/pc.me/pc.c^2;        % classical electron radius [m]
pc.sigma_e  = 8*pi/3*pc.re^2;       % Thomson cross section [m2]
pc.alpha0   = pc.hbar^2/pc.ke/pc.me/pc.e^2;     % Bohr radius [m]
pc.Eh       = pc.alpha^2*pc.c^2*pc.me;          % Hartree energy [J]
pc.Ry       = pc.Eh/2;                          % Rydberg unit of energy [J]
pc.Rinf     = pc.alpha^2*pc.me*pc.c/2/pc.h;     % Rydberg constant [m−1]
pc.Gf       = (pc.h*pc.c)^3;                    % Fermi coupling constant [J−2]
pc.NA       = 6.02214076e23;        % Avogadro constant [mol−1]
pc.R        = pc.NA*pc.kB;          % molar gas constant [J⋅mol−1⋅K−1]
pc.F        = pc.NA*pc.e;           % Faraday constant [C⋅mol−1]
pc.NAh      = pc.NA*pc.h;           % molar Planck constant [J⋅s⋅mol−1]
pc.m_C12    = 1.99264687992e-26;    % atomic mass of carbon-12 [kg]
pc.M_C12    = 11.9999999958E-3;     % molar mass of carbon-12 [kg⋅mol−1]
pc.mu       = 1.66053906660e-27;    % atomic mass constant [kg]
pc.Mu       = 0.99999999965e-3;     % molar mass constant [kg⋅mol−1]
end