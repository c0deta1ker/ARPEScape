function dataStr = load_adress_data(FileName, PathName)
% dataStr = load_adress_data(FileName, PathName)
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
%   dataStr - MATLAB data structure containing all fields below;
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
disp('Loading ADRESS data...')
% - Verifying inputs are characters
FileName = char(FileName);
PathName = char(PathName);

%% 1 - Loading and reading in the HDF5 data file
dataStr = create_arpes_data();
% - If no filename if defined, load in an empty data structure
if isempty(FileName); return;
% - Verify that a .h5 file has been parsed and load the variables
elseif string(FileName(end-2:end)) == ".h5"
    %% 1.1 -- Extracting file information
    FileInfo = dir(char(string(PathName) + string(FileName))); 
    TimeStamp = string(FileInfo.date);
    %% 1.2 -- Extracting all of the data variables
    [Data, Axes, Note]  = ReaderHDF5(char(string(PathName) + string(FileName)));
    Angle               = Axes{1};
    Energy              = (Axes{2})';
    if size(Axes,1)==3; Scan = Axes{3}; else Scan = []; end
    if ndims(Data)==3;  Data = double(permute(Data,[2 1 3])); else Data=double(Data'); end
    % -- Identifying the type of scan that has been performed
    if size(Scan, 1) == 0;          Type = "Eb(k)";         % --- If no Scan parameter is defined, it must be an Eb(k) scan
    elseif range(Scan) == 0;        Type = "Eb(k,i)";       % --- If all Scan parameters are identical, it must be repeated Eb(k) scans
    elseif max(Scan(:)) > 100;      Type = "Eb(kx,kz)";     % --- If the maximum value of the Scan parameter is > 100, it must be Eb(kx,kz)
    else;                           Type = "Eb(kx,ky)";     % --- Else the final type remaining is the Eb(kx,ky) scan
    end
    % -- Extracting the photon energy and mechanical tilt values
    % --- For Eb(k) scans
    if Type == "Eb(k)" || Type == "Eb(k,i)"
        % -- Assigning the photon energy (hv)
        tmpPos = strfind(Note,'hv      =');
        % ---- Vanilla form hv = xxx eV
        for i = 0:1:3
            hv = str2double(Note(tmpPos+9:tmpPos+13-i));
            if ~isnan(hv); break; end
        end
        % ---- If there is a ones(1)* prefix to the hv value
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+20:tmpPos+24-i));
                if ~isnan(hv); break; end
            end
        end
        % ---- If there is a missing space
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+10:tmpPos+14-i));
                if ~isnan(hv); break; end
            end
        end
        % ---- If there is a *ones(1) postfix to the hv value
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+19:tmpPos+23-i));
                if ~isnan(hv); break; end
            end
        end
        if isnan(hv); error('hv in load_data not assigned.'); end
        % -- Assigning the mechanical tilt angle (Tilt)
        tmpPos = strfind(Note, 'Tilt    =');
        Tilt = str2double(Note(tmpPos+9:tmpPos+17));   
        for i = 0:1:7
            Tilt = str2double(Note(tmpPos+9:tmpPos+17-i));   
            if ~isnan(Tilt); break;end
        end
        if isnan(Tilt);error('Tilt in load_data not assigned.'); end
    % --- For Eb(kx,kz) or photon energy dependent scans
    elseif Type == "Eb(kx,kz)"
        % -- Assigning the photon energy (hv) as the scan parameter
        hv = Scan;
        % -- Assigning the tilt angle (Tilt)
        tmpPos = strfind(Note, 'Tilt    =');
        Tilt = str2double(Note(tmpPos+9:tmpPos+17));   
        for i = 0:1:7
            Tilt = str2double(Note(tmpPos+9:tmpPos+17-i));   
            if ~isnan(Tilt); break;end
        end
        if isnan(Tilt);error('Tilt in load_data not assigned.'); end
    elseif Type == "Eb(kx,ky)"
        % -- Assigning the photon energy (hv)
        tmpPos = strfind(Note,'hv      =');
        % -- Vanilla form hv = xxx eV
        for i = 0:1:3
            hv = str2double(Note(tmpPos+9:tmpPos+13-i));
            if ~isnan(hv); break; end
        end
        % -- If there is a ones(1)* prefix to the hv value
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+20:tmpPos+24-i));
                if ~isnan(hv); break; end
            end
        end
        % -- If there is a a missing space
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+10:tmpPos+14-i));
                if ~isnan(hv); break; end
            end
        end
        % -- If there is a *ones(1) postfix to the hv value
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+19:tmpPos+23-i));
                if ~isnan(hv); break; end
            end
        end
        if isnan(hv); error('hv in load_data not assigned.'); end
        % -- Assigning the tilt angle (Tilt) as the scan parameter
        Tilt = Scan;
    end
    %% 1.3 -- Extracting the information variables from the Notes
    % -- Pass energy (ep) and uncertainty (dE) evaluations
    tmpPos = strfind(Note, 'Epass   =');
    for i = 0:1:3
        ep = str2double(Note(tmpPos+9:tmpPos+12-i));
        if ~isnan(ep); break; end
    end
    if isnan(ep); error('ep in load_data not assigned.'); end
    dHv = 75e-3;
    dEnergy = 0.5*ep/1000;
    % -- Theta Manipulator (Theta) evaluation
    tmpPos = strfind(Note, 'Theta   =');
    for i = 0:1:7
        Theta = str2double(Note(tmpPos+9:tmpPos+17-i));
        if ~isnan(Theta); break; end
    end
    if isnan(Theta); error('Theta in load_data not assigned.'); end
    % -- Temperature (Temp) evaluation
    tmpPos = strfind(Note, 'Temp     =');
    for i = 0:1:5
        Temp = str2double(Note(tmpPos+10:tmpPos+14-i));
        if ~isnan(Temp); break; end
    end
    % -- Deflection Angle (Electrostatic Tilt) evaluation
    tmpPos = strfind(Note, 'ADef    =');
    for i = 0:1:5
        ADef = str2double(Note(tmpPos+10:tmpPos+14-i));
        if ~isnan(ADef); break; end
    end
    %% 1.4 - Assigning the data to the MATLAB structure
    % - Assigning file information
    dataStr.FileName    = FileName(1:end-3);
    dataStr.H5file      = char(FileName);
    dataStr.TimeStamp   = TimeStamp;
    dataStr.Type        = Type;
    dataStr.index       = 1:size(Data, 3);
    % - Assigning the meta data
    dataStr.meta.info   = Note;
    dataStr.meta.ep     = ep;
    % - Assinging the main experimental variables
    dataStr.raw_data    = Data;
    dataStr.raw_tht     = Angle;
    dataStr.raw_eb      = Energy;
    dataStr.deb         = dEnergy;
    dataStr.hv          = hv;
    dataStr.dhv         = dHv;
    dataStr.tltM        = Tilt; % Mechanical Tilt
    dataStr.tltE        = ADef; % Electrostatic Tilt
    % - Assinging all other experimental variables
    dataStr.thtM        = Theta;
    dataStr.Temp        = Temp;

%% 2 - Reading in the .mat data that has been processed before
elseif string(FileName(end-3:end)) == ".mat"
    arpes_data  = load(char(string(PathName) + string(FileName)));
    dataStr     = arpes_data.dataStr;
    if isstruct(dataStr) && isfield(dataStr, 'FileName')
        dataStr.FileName = FileName(1:end-4);
    end
end
end