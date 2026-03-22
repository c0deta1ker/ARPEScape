function dataStr_Sum = xsumseries_scans(dataStr)
% dataStr = xcorr_scans(dataStr, scancorr_args)
%   This function cross-correlates the scans over the defined 
%   scan index and converts it from a 3D data set, into a 2D
%   Eb(k) dispersion. The sum is over the scan parameter and works for only
%   Eb(k,i), Eb(kx,ky) or Eb(kx,kz) scans only.
%
%   IN:
%   -   dataStr:            initial data structure of the ARPES data to combine together.
%
%   OUT:
%   -   dataStr:            final combined data structure of ARPES data.

%% - 1 - Initialising input parameters
dataStr_Sum = struct();
ncell       = length(dataStr);
if ncell == 1
    dataStr_Sum.FileName	= dataStr.FileName;
    dataStr_Sum.H5file 	    = dataStr.H5file;
    dataStr_Sum.TimeStamp   = datetime;
    dataStr_Sum.Type        = "Eb(k)";
    dataStr_Sum.index       = dataStr.index;
    % - Assigning the meta data
    dataStr_Sum.meta        = dataStr.meta;
    % - Assinging experimental parameters
    dataStr_Sum.hv          = dataStr.hv;
    dataStr_Sum.tltM        = dataStr.tltM;
    dataStr_Sum.tltE        = dataStr.tltE;
    dataStr_Sum.thtM        = dataStr.thtM;
    dataStr_Sum.Temp   	    = dataStr.Temp;
    % - Assinging the main experimental variables
    dataStr_Sum.raw_data 	= dataStr.raw_data;
    dataStr_Sum.raw_tht  	= dataStr.raw_tht;
    dataStr_Sum.raw_eb   	= dataStr.raw_eb;
    dataStr_Sum.index       = 1;
    dataStr_Sum.raw_data    = sum(dataStr.raw_data, 3);
else
    dataStr_Sum.FileName	= dataStr{1}.FileName;
    dataStr_Sum.H5file 	    = dataStr{1}.H5file;
    dataStr_Sum.TimeStamp   = datetime;
    dataStr_Sum.Type        = "Eb(k)";
    dataStr_Sum.index       = dataStr{1}.index;
    % - Assigning the meta data
    dataStr_Sum.meta        = dataStr{1}.meta;
    % - Assinging experimental parameters
    dataStr_Sum.hv          = dataStr{1}.hv;
    dataStr_Sum.tltM        = dataStr{1}.tltM;
    dataStr_Sum.tltE        = dataStr{1}.tltE;
    dataStr_Sum.thtM        = dataStr{1}.thtM;
    dataStr_Sum.Temp   	    = dataStr{1}.Temp;
    % - Assinging the main experimental variables
    dataStr_Sum.raw_data 	= dataStr{1}.raw_data;
    dataStr_Sum.raw_tht  	= dataStr{1}.raw_tht;
    dataStr_Sum.raw_eb   	= dataStr{1}.raw_eb;
    dataStr_Sum.index       = 1;
    for i = 1:ncell
        raw_data    = cat(3, raw_data, dataStr{i}.raw_data);
    end
    dataStr_Sum.raw_data    = sum(raw_data, 3);
end

end
