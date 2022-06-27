function [imfp, dimfp] = eimfp_optical(ke_dat, element)
% [imfp, dimfp] = eimfp_optical(ke_dat, element)
%   Function that determines the mean and uncertainty of the experimental values
%   of the electron IMFP taken from optical / experimental data. The data
%   used is from the NIST Electron Inelastic-Mean-Free-Path Database [1].
%   Any NaN output values means no experimental data exists for the input
%   kinetic energy, otherwise, the output is the mean and range of the
%   experimental values.
%
%   REFERENCE:
%   [1] NIST Electron Inelastic-Mean-Free-Path Database: http://dx.doi.org/10.18434/T48C78.
%
%   IN:
%   -   ke_dat:     N×1 vector of the input electron kinetic energy (for PES; KE = BE - PHI) [eV]
%   -   element:	char/string of the element; e.g. "H", "He", "Si", "In"...
%
%   OUT:
%   -   imfp:       N×1 column vector of the electron IMFP values [Angstroms]
%   -   dimfp:      N×1 column vector of the IMFP uncertainties [Angstroms]

%% Default parameters (Parameters for Silicon)
if nargin < 2; element = []; end
if isempty(element); element = []; end
%% - 1 - Extracting the IMFP from the optical IMFP database
[ek0, imfp0] = get_eimfpd_props(element);
%% - 2 - Calculate the average and variance of the IMFP values for all datasets
imfp = []; dimfp = [];
% - Filing through all the kinetic energies to be determined
for i = 1:length(ke_dat)
    % -- Extract the value of the closest kinetic energy from the optical data
    ke_i    = ke_dat(i);
    imfp_i  = [];
    for j = 1:size(ek0, 2)
        imfp_i(j)  = interp1(ek0(:,j), imfp0(:,j), ke_i, 'linear');
    end
    % -- Extract the value of the closest kinetic energy from the optical data
    imfp_i(isnan(imfp_i)) = [];
    if isempty(imfp_i)
        imfp(i)     = NaN;
        dimfp(i)    = NaN;
    else
        imfp(i)     = mean(imfp_i);
        dimfp(i)    = 0.5*range(imfp_i);
    end
end
%% Ensuring the imfp is a column vector
if size(imfp, 2) >1; imfp = imfp'; end
if size(dimfp, 2) >1; dimfp = dimfp'; end
end