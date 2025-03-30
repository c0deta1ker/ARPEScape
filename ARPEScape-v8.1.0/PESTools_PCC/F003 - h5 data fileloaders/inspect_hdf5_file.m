function [h5_Info, h5_Data] = inspect_hdf5_file(FileName, PathName)

% -- Open all datasets
h5disp(string(PathName) + string(FileName))


h5_Info = h5info(char(string(PathName) + string(FileName)));
h5_Data         = 0;

% h5read(char(string(PathName) + string(FileName)), '/');



% out_att     = h5readatt(char(string(PathName) + string(FileName)));


end