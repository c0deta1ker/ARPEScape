function pesStr = create_pes_data()
% pesStr = create_pes_data()
%   This function generates a template data structure used for
%   analysing PES data. The fields can be filled in, after which the user
%   can then use the data structure within all the processing / fitting
%   tools.
pesStr          = struct();
pesStr.FileName	= '';
pesStr.H5file 	= '';
pesStr.Type     = "PES";
pesStr.hv       = [];
pesStr.xdat     = [];
pesStr.ydat   	= [];
end