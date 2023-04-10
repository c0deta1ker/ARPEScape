function [Ec0, Ev0] = semi_si_nd2ebands(nd)
% [Ec0, Ev0] = semi_si_nd2ebands(nd)
%   Function that determines the CBM / VBM energy for a given dopant
%   concentration. The data is taken from [1] using Phosphorous as the
%   dopant. The band gap narrowing has been implemented and the Fermi-level
%   is located at 0 eV (i.e. the output values of the CBM / VBM are given
%   relative to the Fermi-level).
%
%   REFERENCE:
%   [1] PV Lighthouse Solar Cell Simulator: https://www.pvlighthouse.com.au/bandgap
%
%   IN:
%   -   nd:         N×1 column vector of the dopant concentration [cm-3]
%
%   OUT:
%   -   Ec0:        N×1 column vector of CBM relative to EF [eV]
%   -   Ev0:        N×1 column vector of VBM relative to EF [eV]

%% Default parameters
if nargin < 1; nd = 0; end
if isempty(nd); nd = 0; end
%% - 1 - Extracting the energy bands relative to EF vs dopant concentration
% - Importing the relevant data file
raw_dat = importdata('si_ebands_vs_nd.txt');
nd_dat  = raw_dat.data(:,1);    % dopant concentration
Ec0_dat = raw_dat.data(:,2);    % CBM energy (eV)
Ev0_dat = raw_dat.data(:,2);    % VBM energy (eV)
% - Linear interpolation of the data
Ec0_F   = fit(nd_dat, Ec0_dat, 'linearinterp');
Ev0_F   = fit(nd_dat, Ev0_dat, 'linearinterp');
% - Extracting the CBM and VBM values in eV
Ec0     = Ec0_F(nd);
Ev0     = Ev0_F(nd);
end