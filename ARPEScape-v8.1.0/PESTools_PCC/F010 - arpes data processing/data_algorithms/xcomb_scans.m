function dataStr = xcomb_scans(varargin)
% dataStr = xcorr_scans(dataStr, scancorr_args)
%   This function cross-correlates the scans over the defined 
%   scan index and converts it from a 3D data set, into a 2D
%   Eb(k) dispersion. The sum is over the scan parameter and works for only
%   Eb(k,i), Eb(kx,ky) or Eb(kx,kz) scans only.
%
%   IN:
%   -   varargin:           arguments: each data structure of the ARPES data to combine together.
%
%   OUT:
%   -   dataStr:            fianl combined data structure of ARPES data.

%% - 1 - Initialising input parameters
if ~isempty(varargin)
    dataStr.FileName	= varargin{1}.FileName;
    dataStr.H5file 	    = varargin{1}.H5file;
    dataStr.TimeStamp   = datetime;
    dataStr.Type        = "Eb(k,i)";
    dataStr.index       = varargin{1}.index;
    % - Assigning the meta data
    dataStr.meta        = varargin{1}.meta;
    % - Assinging experimental parameters
    dataStr.hv          = varargin{1}.hv;
    dataStr.tltM        = varargin{1}.tltM;
    dataStr.tltE        = varargin{1}.tltE;
    dataStr.thtM        = varargin{1}.thtM;
    dataStr.Temp   	    = varargin{1}.Temp;
    % - Assinging the main experimental variables
    dataStr.raw_data 	= varargin{1}.raw_data;
    dataStr.raw_tht  	= varargin{1}.raw_tht;
    dataStr.raw_eb   	= varargin{1}.raw_eb;
    if length(varargin) > 1
        for i = 2:length(varargin)
            dataStr.index       = horzcat(dataStr.index, dataStr.index(end)+1 : 1 : dataStr.index(end) + length(varargin{i}.index));
            dataStr.raw_data    = cat(3, dataStr.raw_data, varargin{i}.raw_data);
        end
    end
end
end
