function nd = semi_si_Ev02nd(Ev0)
% nd = semi_si_Ev02nd(Ev0)
%   Function that determines the dopant concentration based on the VBM
%   position relative to the Fermi-level.The data is taken from [1] using 
%   Phosphorous as the dopant.
%
%   REFERENCE:
%   [1] PV Lighthouse Solar Cell Simulator: https://www.pvlighthouse.com.au/bandgap
%
%   IN:
%   -   Ev0:        N×1 column vector of VBM relative to EF [eV] 
%
%   OUT:
%   -   nd:         N×1 column vector of the dopant concentration [cm-3]

%% Default parameters
if nargin < 1; Ev0 = 0; end
if isempty(Ev0); Ev0 = 0; end
%% - 1 - Extracting the energy bands relative to EF vs dopant concentration
% - Importing the relevant data file
raw_dat = importdata('si_ebands_vs_nd.txt');
nd_dat  = raw_dat.data(:,1);    % dopant concentration
Ec0_dat = raw_dat.data(:,2);    % CBM energy (eV)
% - Linear interpolation of the data
nd_F   = fit(Ec0_dat, nd_dat,'linearinterp');
% - Extracting the dopant concentration
nd     = nd_F(Ev0);
end