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
%
%   IN:
%   -   FileName:           char of the input .h5 or .mat file-name.
%   -   PathName:           char of the input .h5 or .mat directory path.
%
%   OUT:
%   dataStr - MATLAB data structure containing all fields below;
%	 .(FileName):       string of the current filename of the data file.
%	 .(H5file):         string of the raw .H5 filename of the data file.
%	 .(TimeStamp):      time-stamp of when the file was created.
%	 .(Type):           string "Eb(k)", "Eb(k,i)", "Eb(kx,ky)" or "Eb(kx,kz)".
% 	 .(index):          [1xN] vector of the total number of scans.
%	 .(meta):           structure that contains all meta information.
%	    .(meta.info):       [1xN] char array of scan information.
%	    .(meta.Epass):      scalar of the Pass energy.
%	    .(meta.Pol):        string of the photon polarisation used.
%	    .(meta.Slit):       scalar of the beamline exit slit size.
%	    .(meta.Mode):       string of the Analyzer angular mode used.
%	    .(meta.X):          scalar of the manipulator X position.
%	    .(meta.Y):          scalar of the manipulator Y position.
%	    .(meta.Z):          scalar of the manipulator Z position.
%	    .(meta.dhv):        scalar estimate of the beamline resolution.
%	    .(meta.deb):        scalar of the Analyzer resolution.
%	 .(hv):              scalar or [1×nScan] row vector (if Scan parameter for 3Data data). Photon energy [eV].
%	 .(tltM):            scalar or [1×nScan] row vector (if Scan parameter for 3Data data). Mechanical Tilt [degrees]
%	 .(tltE):            scalar or [1×nScan] row vector (if Scan parameter for 3Data data). Electrostatic Tilt [degrees].
%	 .(thtM):            scalar of the manipulator Theta angle.
%	 .(Temp):            scalar of the sample temperature.
%	 .(raw_data):        2Data [nEnergy x nAngle] or 3Data [nEnergy x nAngle xnScan] data array.
%	 .(raw_tht):         [1×nAngle] row vector of angles.
% 	 .(raw_eb):          [nEnergyx1] column vector of energy.


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
    FileInfo    = dir(char(string(PathName) + string(FileName))); 
    TimeStamp   = string(FileInfo.date);
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
        % ---- Vanilla form hv = xxxx eV
        for i = 0:1:3
            hv = str2double(Note(tmpPos+9:tmpPos+13-i));
            if ~isnan(hv); break; end
        end
        % ---- If there is a missing space
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+10:tmpPos+14-i));
                if ~isnan(hv); break; end
            end
        end
        % ---- If there is a ones(x)* prefix to the hv value
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+18:tmpPos+24-i));
                if ~isnan(hv); break; end
            end
        end
        % ---- If there is a ones(xx)* prefix  to the hv value
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+19:tmpPos+23-i));
                if ~isnan(hv); break; end
            end
        end
        if isnan(hv); hv = []; disp('hv in load_adress_data not assigned.'); end
        % -- Assigning the mechanical tilt angle (Tilt)
        tmpPos = strfind(Note, 'Tilt    =');
        Tilt = str2double(Note(tmpPos+9:tmpPos+17));   
        for i = 0:1:7
            Tilt = str2double(Note(tmpPos+9:tmpPos+17-i));   
            if ~isnan(Tilt); break;end
        end
        if isnan(Tilt); Tilt = []; disp('Tilt in load_adress_data not assigned.'); end
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
        if isnan(Tilt); Tilt = []; disp('Tilt in load_adress_data not assigned.'); end
    elseif Type == "Eb(kx,ky)"
        % -- Assigning the photon energy (hv)
        tmpPos = strfind(Note,'hv      =');
        % ---- Vanilla form hv = xxxx eV
        for i = 0:1:3
            hv = str2double(Note(tmpPos+9:tmpPos+13-i));
            if ~isnan(hv); break; end
        end
        % ---- If there is a missing space
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+10:tmpPos+14-i));
                if ~isnan(hv); break; end
            end
        end
        % ---- If there is a ones(x)* prefix to the hv value
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+18:tmpPos+24-i));
                if ~isnan(hv); break; end
            end
        end
        % ---- If there is a ones(xx)* prefix  to the hv value
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+19:tmpPos+23-i));
                if ~isnan(hv); break; end
            end
        end
        if isnan(hv); hv = []; disp('hv in load_adress_data not assigned.'); end
        % -- Assigning the tilt angle (Tilt) as the scan parameter
        Tilt = Scan;
    end
    %% 1.3 -- Extracting the information variables from the Notes
    % -- Pass energy (ep) and uncertainty (dE) evaluations
    tmpPos = strfind(Note, 'Epass   =');
    for i = 0:1:3
        Epass = str2double(Note(tmpPos+9:tmpPos+12-i));
        if ~isnan(Epass); break; end
    end
    if isnan(Epass); Epass = []; disp('Epass in load_adress_data not assigned.'); end
    dHv = 75e-3;
    dEnergy = 0.5*Epass/1000;
    % -- Theta Manipulator (Theta) evaluation
    tmpPos = strfind(Note, 'Theta   =');
    for i = 0:1:7
        Theta = str2double(Note(tmpPos+9:tmpPos+17-i));
        if ~isnan(Theta); break; end
    end
    if isnan(Theta); Theta = []; disp('Theta in load_adress_data not assigned.'); end
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
    % -- Polarisation evaluation
    tmpPos = strfind(Note, 'Pol     =');
    for i = 0:1:5
        Pol = string(Note(tmpPos+10:tmpPos+12-i));
        if contains(Pol,"LV"); Pol = "LV (p-pol)"; break;
        elseif contains(Pol,"LH"); Pol = "LH (s-pol)"; break;
        elseif contains(Pol,"C+"); Pol = "C+"; break;
        elseif contains(Pol,"C-"); Pol = "C-"; break;
        else; Pol = "";
        end
    end
    % -- Slit evaluation
    tmpPos = strfind(Note, 'Slit    = ');
    for i = 0:1:5
        Slit = str2double(Note(tmpPos+10:tmpPos+14-i));
        if ~isnan(Slit); break; end
    end
    % -- Mode evaluation
    tmpPos = strfind(Note, 'Mode    = ');
    for i = 0:1:5
        Mode = string(Note(tmpPos+10:tmpPos+13-i));
        if contains(Mode,"MAD"); Mode = "MAD"; break;
        elseif contains(Mode,"MAM"); Mode = "MAM"; break;
        elseif contains(Mode,"LAD"); Mode = "LAD"; break;
        elseif contains(Mode,"WAM"); Mode = "WAM"; break;
        else; Mode = "";
        end
    end
    % -- Manipulator Position evaluation
    % --- X
    tmpPos = strfind(Note, 'X       = ');
    for i = 0:1:8
        X = str2double(Note(tmpPos+10:tmpPos+20-i));
        if ~isnan(X); break; end
    end
    % --- Z
    tmpPos = strfind(Note, 'Z       = ');
    for i = 0:1:8
        Z = str2double(Note(tmpPos+10:tmpPos+18-i));
        if ~isnan(Z); break; end
    end
    % --- Y
    tmpPos = strfind(Note, 'Y        = ');
    for i = 0:1:6
        Y = str2double(Note(tmpPos+11:end-i));
        if ~isnan(Y); break; end
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
    dataStr.meta.Epass  = Epass;
    dataStr.meta.Pol    = Pol;
    dataStr.meta.Slit   = Slit;
    dataStr.meta.Mode   = Mode;
    dataStr.meta.X      = X;
    dataStr.meta.Y      = Y;
    dataStr.meta.Z      = Z;
    dataStr.meta.dhv    = dHv;
    dataStr.meta.deb    = dEnergy;
    % - Assinging experimental parameters
    dataStr.hv          = hv;
    dataStr.tltM        = Tilt;     % Mechanical Tilt
    dataStr.tltE        = ADef;     % Electrostatic Tilt
    dataStr.thtM        = Theta;
    dataStr.Temp        = Temp;
    % - Assinging the main experimental variables
    dataStr.raw_data    = Data;
    dataStr.raw_tht     = Angle;
    dataStr.raw_eb      = Energy;

%% 2 - Reading in the .mat data that has been processed before
elseif string(FileName(end-3:end)) == ".mat"
    arpes_data  = load(char(string(PathName) + string(FileName)));
    dataStr     = arpes_data.dataStr;
    if isstruct(dataStr) && isfield(dataStr, 'FileName')
        dataStr.FileName = FileName(1:end-4);
    end
end
end