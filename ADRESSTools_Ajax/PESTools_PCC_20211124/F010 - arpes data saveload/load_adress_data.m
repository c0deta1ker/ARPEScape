function dataStr = load_adress_data(FileName, PathName, Reset)
% dataStr = load_adress_data(FileName, PathName, Reset)
%   This function loads in the HDF5 data files of ARPES data from 
%   the ADRESS beamline at the SLS. The output is a MATLAB
%   data-structure that yields all the data and information. This function 
%   should be used to load in a single ADRESS data-file that is unprocessed
%   or has been previously processed using PESTools.
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
% -- after processing (1) Eb alignment
%	 .(eb_ref_file):     string of the reference file-name used (if used).
%	 .(eb_shifts):       cell-array of the binding energy shifts.
%	 .(eb):              aligned 2D or 3D array of energy.
%	 .(tht):             aligned 2D or 3D array of theta.
% -- after processing (2) intensity normalisation
%	 .(data):            2D or 3D array of the normalised data.
% -- after processing (3) k conversions
%	 .(surfNormX):       double or vector of surface normal vector.
%	 .(kx):              2D or 3D array of kx from the Theta angle.
%	 .(ky):              1D, 2D or 3D array of kx from the Tilt angle.
%	 .(kz):              1D, 2D or 3D array of kx from the Photon Energy.

%% Default parameters
if nargin < 2; PathName = ''; Reset = 0; end
if nargin < 3; Reset = 0; end
if isempty(PathName);   PathName = '';  end
if isempty(Reset);      Reset = 0;  end
disp('Loading ARPES data...')
% wbar = waitbar(0, 'Loading in ARPES data...', 'Name', 'load_adress_data');
% wbar.Children.Title.Interpreter = 'none';
% - Verifying input are characters
FileName = char(FileName);
PathName = char(PathName);

%% 1 - Reading in the HDataF5 data
% - Verify that a .h5 file has been parsed and load the variables
if string(FileName(end-2:end)) == ".h5"
%     waitbar(0.1, wbar, sprintf('Loading %s...', FileName), 'Name', 'load_adress_data');
    % -- Defining the MATLAB data-structure
    dataStr             = struct();
    dataStr.FileName    = FileName(1:end-3);
    dataStr.H5file      = char(FileName);
    % -- Extracting all of the data variables
    [Data, Axes, Note]  = ReaderHDF5(char(string(PathName) + string(FileName)));
    Angle               = Axes{1};
    Energy              = (Axes{2})';
    if size(Axes,1)==3; Scan = Axes{3}; else Scan = []; end
    if ndims(Data)==3;  Data = double(permute(Data,[2 1 3])); else Data=double(Data'); end
    % -- Identifying the type of scan that has been performed
    % --- If no Scan parameter is initially defined, it is an Eb(k) scan
    if size(Scan, 1) == 0
        dataStr.Type = "Eb(k)";
        % -- Assigning the photon energy (hv)
        tmpPos = strfind(Note,'hv      =');
        for i = 0:1:3
            hv = str2double(Note(tmpPos+9:tmpPos+13-i));
            if ~isnan(hv); break; end
        end
        % -- If there is a ones(1)* prefix to the hv value
        if isnan(hv)
            for i = 0:1:3
                hv = str2double(Note(tmpPos+18:tmpPos+22-i));
                if ~isnan(hv); break; end
            end
        end
        if isnan(hv); error('hv in load_data not assigned.'); end
        % -- Assigning the tilt angle (Tilt)
        tmpPos = strfind(Note, 'Tilt    =');
        Tilt = str2double(Note(tmpPos+9:tmpPos+17));   
        for i = 0:1:7
            Tilt = str2double(Note(tmpPos+9:tmpPos+17-i));   
            if ~isnan(Tilt); break;end
        end
        if isnan(Tilt);error('Tilt in load_data not assigned.'); end
    % --- If the first and last Scan parameter are identical, multiple Eb(k) scans have been performed
    elseif Scan(1) - Scan(end) == 0
        dataStr.Type = "Eb(k)";
        % -- Assigning the photon energy (hv)
        tmpPos = strfind(Note,'hv      =');
        for i = 0:1:3
            hv = str2double(Note(tmpPos+9:tmpPos+13-i));
            if ~isnan(hv); break; end
        end
        if isnan(hv); error('hv in load_data not assigned.'); end
        % -- Assigning the tilt angle (Tilt)
        tmpPos = strfind(Note, 'Tilt    =');
        Tilt = str2double(Note(tmpPos+9:tmpPos+17));   
        for i = 0:1:7
            Tilt = str2double(Note(tmpPos+9:tmpPos+17-i));   
            if ~isnan(Tilt); break;end
        end
        if isnan(Tilt);error('Tilt in load_data not assigned.'); end
    % --- If the maximum value of the Scan parameter is > 50, it must be Eb(kx,kzs)
    elseif max(Scan(:)) > 50
        dataStr.Type = "Eb(kx,kz)";
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
    % --- Else the final type remaining is the Eb(kx,ky) scan
    else
        dataStr.Type = "Eb(kx,ky)";
        % -- Assigning the photon energy (hv)
        tmpPos = strfind(Note,'hv      =');
        for i = 0:1:3
            hv = str2double(Note(tmpPos+9:tmpPos+13-i));
            if ~isnan(hv); break; end
        end
        if isnan(hv); error('hv in load_data not assigned.'); end
        % -- Assigning the tilt angle (Tilt) as the scan parameter
        Tilt = Scan;
    end
    %% 1.2 - Extracting the information variables from the Notes
%     waitbar(0.25, wbar, 'Extracting information variables...', 'Name', 'load_adress_data');
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
    %% 1.3 - For an Eb(k) scan with multiple scans, cross-correlate them
    if dataStr.Type == "Eb(k)" && size(Data, 3) > 1
        dataStr.Type = "Eb(k,i)";
        dataStr.index = 1:size(Data, 3);
%         waitbar(0.50, wbar, 'Cross-correlating repeated scans...', 'Name', 'load_adress_data');
%         [xc_Data, ~]            = SumScanXC(Energy, Data, 0.5); 
%         xc_Data(isnan(xc_Data)) = 0;
%         Data                    = xc_Data;
    end
    %% 1.4 - Assigning the data to the MATLAB structure
%     waitbar(0.75, wbar, 'Assigning ARPES data to MATLAB structure...', 'Name', 'load_adress_data');
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
    dataStr.tltM        = Tilt;
    % - Assinging all other experimental variables
    dataStr.thtM        = Theta;
    dataStr.Temp        = Temp;

%% 2 - Reading in the .mat data that has been processed
elseif string(FileName(end-3:end)) == ".mat"
%     waitbar(0.5, wbar, sprintf('Loading %s...', FileName), 'Name', 'load_adress_data');
    arpes_data  = load(char(string(PathName) + string(FileName)));
    % Loading data to a MATLAB structure
    dataStr     = arpes_data.dataStr;
    if isstruct(dataStr) && isfield(dataStr, 'FileName')
        dataStr.FileName = FileName(1:end-4);
    end
end

%% 3 - Removing all other processing that may have been performed in other apps
if Reset == 1
    % -- Kf analysis
    if isfield(dataStr, 'kf');      dataStr = rmfield(dataStr, 'kf'); end
    % -- Isoe analysis
    if isfield(dataStr, 'isoe');    dataStr = rmfield(dataStr, 'isoe'); end
    if isfield(dataStr, 'bz');      dataStr = rmfield(dataStr, 'bz'); end
    % -- State-fitting analysis
    if isfield(dataStr, 'fits');    dataStr = rmfield(dataStr, 'fits'); end
end

%% Close wait-bar
% close(wbar);
end