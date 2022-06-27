function arpesStr = create_arpes_data()
% arpesStr = create_arpes_data()
%   This function generates a template data structure used for
%   analysing ARPES data. The fields can be filled in, after which the user
%   can then use the data structure within all the processing / fitting
%   tools.
arpesStr            = struct();
arpesStr.FileName	= '';
arpesStr.H5file 	= '';
arpesStr.Type       = "Eb(k)";
arpesStr.index      = 1;
arpesStr.meta       = struct();
arpesStr.raw_data 	= [];
arpesStr.raw_tht  	= [];
arpesStr.raw_eb   	= [];
arpesStr.deb      	= [];
arpesStr.hv         = [];
arpesStr.dhv        = [];
arpesStr.tltM     	= [];
arpesStr.thtM       = [];
arpesStr.Temp   	= [];
end