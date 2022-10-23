function pearlStr = load_pearl_data(FileName, PathName)
% pearlStr = load_adress_data(FileName, PathName)
%   This function loads in the HDF5 data files of ARPES data from 
%   the ADRESS beamline at the SLS. The output is a MATLAB
%   data-structure that yields all the data and information. This function 
%   should be used to load in a single ADRESS data-file that is unprocessed
%   or has been previously processed using PESTools. Last updated in
%   October 2022, and includes electrostatic tilt extraction.
%
%   REQ. FUNCTIONS:
%   -   [data, axes, note]  = ReaderHDF5(fname);
%   -   [DataXC,OffsE]      = SumScanXC(Energy,Data,maxLagE[,WinE]);
%
%   IN:
%   -   FileName:           char of the input .h5 or .mat file-name.
%   -   PathName:           char of the input .h5 or .mat directory path.
%   -   Reset:              if equal to 1, then reset all previous processing upon loading, otherwise leave alone.
%
%   OUT:
%   pearlStr - MATLAB data structure containing all fields below;
%	 .(FileName):        string of the current filename of the data file.
%	 .(H5file):          string of the raw .H5 filename of the data file.
%	 .(Type):            string "Eb(k)", "Eb(k,i)", "Eb(kx,ky)" or "Eb(kx,kz)".
%	 .(meta.info):       [1xN] char array of scan information.
%	 .(meta.ep):         scalar of the Pass energy.
%	 .(raw_data):        2Data [nEnergy x nAngle] or 3Data [nEnergy x nAngle xnScan] data array.
%	 .(raw_tht):         [1×nAngle] row vector of angles.
% 	 .(raw_eb):          [nEnergyx1] column vector of energy.
%	 .(deb):             scalar of the Analyser resolution.
%	 .(hv):              scalar or [1×nScan] row vector (if Scan parameter for 3Data data).
%	 .(dhv):             scalar of the Beamline resolution.
%	 .(tltM):            scalar or [1×nScan] row vector (if Scan parameter for 3Data data).
%	 .(thtM):            scalar of the manipular Theta angle.
%	 .(Temp):            scalar of the sample temperature.

%% Default parameters
if nargin < 1; FileName = ''; end
if nargin < 2; PathName = ''; end
if isempty(FileName);   FileName = '';  end
if isempty(PathName);   PathName = '';  end
% - Display text
disp('Loading PEARL data...')
% - Verifying inputs are characters
FileName = char(FileName);
PathName = char(PathName);

%% 1 - Loading and reading in the HDF5 data file
% - If no filename if defined, load in an empty data structure
if isempty(FileName); return;
% - Verify that a .h5 file has been parsed and load the variables
elseif string(FileName(end-2:end)) == ".h5"
    %% 1.1 -- Extracting file information
    FileInfo = dir(char(string(PathName) + string(FileName))); 
    TimeStamp = string(FileInfo.date);
    %% 1.2 -- Extracting all of the miscellaneous variables
    pearlStr = struct();
    full_h5_file_name = char(string(PathName) + string(FileName));
    % (1) Finding all the HDF5 Group Names
    Type = "Spectrum";
    group_names = string({h5info(full_h5_file_name).Groups.Name});
    if group_names(1) == '/__DATA_TYPES__'; group_names(1) = []; end
    % -- Filing through each group and extracting all the data
    for i = 1:length(group_names)
        % Validity check on group name
        % ---- Ensuring the group name does not include prefix '/'
        group_name = group_names(i); group_name = char(group_name); group_name = group_name(2:end);
        % ---- Removing all spaces from the group name
        group_name = strrep(group_name," ","_");
        % Extracting element names
        element_names = string({h5info(full_h5_file_name, char(group_names(i)+"/")).Datasets.Name});
        % --- Filing through each element in the group
        for j = 1:length(element_names)
            % Validity check on element name
            element_name = element_names(j);
            % ---- Extracting the PEARL dataset
            pearlStr.(group_name).(element_name) = h5read(full_h5_file_name, char(group_names(i)+"/"+element_names(j)));
        end
    end
    %% 1.3 -- Extracting all of the data variables - 'scan 1'
    if contains(group_names, "/scan 1")
        % For custom 'MultiRegionScan' files
        Type = "MultiRegionScan";
        group_names = string({h5info(full_h5_file_name, char("/scan 1")).Groups.Name});
        % -- Filing through each group and extracting all the data
        for i = 1:length(group_names)
            % ---- Ensuring the group name does not include prefix '/'
            group_name = group_names(i); group_name = char(group_name); group_name = group_name(9:end);
            % Extracting element names
            element_names = string({h5info(full_h5_file_name, char(group_names(i)+"/")).Datasets.Name});
            % --- Filing through each element in the group
            for j = 1:length(element_names)
                pearlStr.scan_1.(group_name).(element_names(j)) = h5read(full_h5_file_name, char(group_names(i)+"/"+element_names(j)));
            end
        end
    end

%% 2 - Reading in the .mat data that has been processed before
elseif string(FileName(end-3:end)) == ".mat"
    arpes_data  = load(char(string(PathName) + string(FileName)));
    pearlStr     = arpes_data.pearlStr;
    if isstruct(pearlStr) && isfield(pearlStr, 'FileName')
        pearlStr.FileName = FileName(1:end-4);
    end
end
end