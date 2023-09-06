function arpesStr = create_arpes_data()
% arpesStr = create_arpes_data()
%   This function generates a template data structure used for
%   analysing ARPES data. The fields can be filled in, after which the user
%   can then use the data structure within all the processing / fitting
%   tools.
% - Assigning file information
arpesStr            = struct();
arpesStr.FileName	= '';
arpesStr.H5file 	= '';
arpesStr.TimeStamp   = datetime;
arpesStr.Type       = "Eb(k)";
arpesStr.index      = 1;
% - Assigning the meta data
arpesStr.meta       = struct();
% - Assinging experimental parameters
arpesStr.hv         = [];
arpesStr.tltM     	= [];
arpesStr.tltE     	= [];
arpesStr.thtM       = [];
arpesStr.Temp   	= [];
% - Assinging the main experimental variables
arpesStr.raw_data 	= [];
arpesStr.raw_tht  	= [];
arpesStr.raw_eb   	= [];
end