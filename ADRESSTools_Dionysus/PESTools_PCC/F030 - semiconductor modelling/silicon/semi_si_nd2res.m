function resistivity = semi_si_nd2res(nd)
% resistivity = semi_si_nd2res(nd)
%   Function that determines the resistivity of a silicon sample,
%   determined by the dopant concentration only.
%
%   REFERENCE:
%   [1] PV Lighthouse Solar Cell Simulator: https://www.pvlighthouse.com.au/resistivity
%
%   IN:
%   -   Ec0:        N×1 column vector of CBM relative to EF [eV] 
%
%   OUT:
%   -   nd:         N×1 column vector of the dopant concentration [cm-3]

%% Default parameters
if nargin < 1; nd = 1e10; end
if isempty(nd); nd = 1e10; end
%% - 1 - Extracting the energy bands relative to EF vs dopant concentration
% - Importing the relevant data file
raw_dat             = importdata('si_nd_vs_res.txt');
nd_dat              = raw_dat.data(:,1);            % dopant concentration
resistivity_dat     = raw_dat.data(:,2);            % resistivity
% - Linear interpolation of the data
resistivity_F   = fit(nd_dat, resistivity_dat, 'linearinterp');
% - Extracting the dopant concentration
resistivity     = resistivity_F(nd);
end