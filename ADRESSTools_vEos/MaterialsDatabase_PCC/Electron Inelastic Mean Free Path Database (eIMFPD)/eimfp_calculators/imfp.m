function imfp = imfp(formalism, ke_dat, args)
% imfp = imfp(formalism, ke_dat, args)
%   Function that calculates the electron inelastic mean free path (IMFP).
%   The TPP2M, S1, S2 and Universal formalisms are available. The user can
%   define a scalar or vector of kinetic energies for the input.
%   If the args is a string of an element / material, it will look it up 
%   the relevant materials parameters in the Material Properties Database 
%   (MPD) ('MPD_PCC.mat'). Otherwise, the args a be manually inserted as
%   discussed below.
%
%   IN:
%   -   formalism:  string for imfp calculator formalism. Default:"TPP2M" ["TPP2M","S1","S2","Universal"]
%   -   ke_dat:  	N×1 column vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%   -   args:       string or vector of scalars defining the arguments of the chosen formalism;
%                       (1) string of element / material whose parameters are found in materials database; e.g. "Si", "SiO2", "Al2O3"...
%                       (2) vector with manual entry of material parameters:
%                            -> TPP2M:  4x1     [density(g/cc),atomicweight(amu),egap(eV),valency(valence electrons per atom)]
%                            -> S1:     4x1     [density(g/cc),atomicweight(amu),egap(eV),Z(atomic mass number or average for compounds)]
%                            -> S2:     1x1     [Z(atomic mass number or average for compounds)]
%                            -> Universal:      no args, material independent.
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]

%% Default parameters
% - Default formalism
if nargin < 3; args = []; end
if nargin < 1; formalism = "TPP2M"; end
if isempty(formalism); formalism = "TPP2M"; end
if isempty(args); args = []; end
% - Validity check on the inputs
if ischar(formalism); formalism = string(formalism); end
if ischar(args); args = string(args); end

%% - 1 - Determination of the IMFP
% If the kinetic energy is negative, assume it is zero
ke_dat(ke_dat<0) = 0;
% -- TPP2M formalism
if strcmpi(formalism, "TPP2M") || strcmpi(formalism, "TPP-2M")
    if isstring(args)
        material = args;
        imfp = eimfp_tpp2m_mpd(ke_dat, material);
    else
        rho = args(1); Nv=args(4); M = args(2); Egap = args(3);
        imfp = eimfp_tpp2m(ke_dat, rho, Nv, M, Egap);
    end
% -- S1 formalism
elseif strcmpi(formalism, "S1")
    if isstring(args)
        material = args;
        imfp = eimfp_S1_mpd(ke_dat, material);
    else
        rho = args(1); Z=args(4); M = args(2); Egap = args(3);
        imfp = eimfp_S1(ke_dat, rho, M, Egap, Z);
    end
% -- S2 formalism
elseif strcmpi(formalism, "S2")
    if isstring(args)
        material = args;
        imfp = eimfp_S2_avg_mpd(ke_dat, material);
    else
        Z=args(1);
        imfp = eimfp_S2_avg(ke_dat, Z);
    end
% -- Universal formalism
elseif strcmpi(formalism, "Universal")
    imfp = eimfp_universal(ke_dat);
end
%% Ensuring the imfp is a column vector
if size(imfp, 2) >1; imfp = imfp'; end
end