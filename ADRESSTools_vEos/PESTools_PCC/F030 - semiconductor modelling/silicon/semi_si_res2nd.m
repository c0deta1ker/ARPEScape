function nd = semi_si_res2nd(resistivity)
% nd = semi_si_res2nd(resistivity)
%   Function that determines the dopant concentration of a silicon sample,
%   determined by the resistivity only. The data is taken from [1] using 
%   Phosphorous as the dopant.
%
%   REFERENCE:
%   [1] PV Lighthouse Solar Cell Simulator: https://www.pvlighthouse.com.au/resistivity
%
%   IN:
%   -   resistivity:        N×1 column vector of the resistivity [Ohm cm] 
%
%   OUT:
%   -   nd:                 N×1 column vector of the dopant concentration [cm-3]

%% Default parameters
if nargin < 1; resistivity = 0; end
if isempty(resistivity); resistivity = 0; end
%% - 1 - Extracting the energy bands relative to EF vs dopant concentration
% - Importing the relevant data file
raw_dat = importdata('si_nd_vs_res.txt');
nd_dat  = raw_dat.data(:,1);            % dopant concentration
resistivity_dat = raw_dat.data(:,2);    % resistivity
% - Linear interpolation of the data
nd_F   = fit(resistivity_dat, nd_dat,'linearinterp');
% - Extracting the dopant concentration
nd     = nd_F(resistivity);
end