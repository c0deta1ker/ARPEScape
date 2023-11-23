function dataStr = load_SES_zip(FileName,PathName)
%% Default parameters
if nargin < 1; PathName = ''; end
% - Display text
disp('Loading SES zip data...')

% - Verifying inputs are characters
FileName = char(FileName);
PathName = char(PathName);
dataStr = create_arpes_data();

% Unzip zip file
full_file_name = char(string(PathName) + string(FileName));
tempDir = fullfile(PathName, 'temp'); % Creates a full path to the 'temp' directory
z = unzip(full_file_name,tempDir);

%% Get file names

% Get region names
for i = 1:numel(z)
    % Individual file names
    [~, name, ext] = fileparts(z{i});
    fullName = [name ext]; % restructure file names

    % check if it starts with 'Spectrum_' and ends with '.bin'
    if startsWith(fullName, 'Spectrum_') && endsWith(fullName, '.bin')
        % delete 'Spectrum_' prefix '.bin' 
        regionName = erase(fullName, 'Spectrum_');
        regionName = erase(regionName, '.bin');
    end
end

% file names
filename_spectrum = ['Spectrum_' regionName '.ini'];
spectrumFile = fullfile(tempDir, filename_spectrum);

filename = [regionName '.ini'];
regionFile = fullfile(tempDir, filename);

filename_bin = ['Spectrum_' regionName '.bin'];
spectrumFile_bin = fullfile(tempDir, filename_bin);

%% Load .ini files
% Open the spectrumFile file
fileID = fopen(spectrumFile, 'r');

% Check if the file is successfully opened
if fileID == -1
    error('File cannot be opened');
end

% Read the file line by line
while ~feof(fileID)
    line = fgetl(fileID);

    %---- Energy axis ----
    if startsWith(line, "widthoffset=")
        LowEnergy = str2double(strrep(strrep(line, "widthoffset=", ""), ",", "."));
    end
    if startsWith(line, "width=")
        numEnergyPoints = str2double(strrep(strrep(line, "width=", ""), ",", "."));
    end
    if startsWith(line, "widthdelta=")
        EnergyStep = str2double(strrep(strrep(line, "widthdelta=", ""), ",", "."));
    end

    %---- Analyzer slit axis ----
    if startsWith(line, "heightoffset=")
        LowAnalyzerAngle = str2double(strrep(strrep(line, "heightoffset=", ""), ",", "."));
    end
    if startsWith(line, "height=")
        numAnalyzerPoints = str2double(strrep(strrep(line, "height=", ""), ",", "."));
    end
    if startsWith(line, "heightdelta=")
        AnalyzerStep = str2double(strrep(strrep(line, "heightdelta=", ""), ",", "."));
    end

    %---- Deflector axis ----
    if startsWith(line, "depthoffset=")
        LowDeflectorAngle = str2double(strrep(strrep(line, "depthoffset=", ""), ",", "."));
    end
    if startsWith(line, "depth=")
        numDeflectorPoints = str2double(strrep(strrep(line, "depth=", ""), ",", "."));
    end
    if startsWith(line, "depthdelta=")
        DeflectorStep = str2double(strrep(strrep(line, "depthdelta=", ""), ",", "."));
    end
end

% Close the file
fclose(fileID);

% Open the regionFile file
fileID = fopen(regionFile, 'r');
% Check if the file is successfully opened
if fileID == -1
    error('File cannot be opened');
end

% Read the file line by line
while ~feof(fileID)
    line = fgetl(fileID);

    % Process each line based on its starting text
    % Extracting metadata from the file based on line prefixes
    if startsWith(line, "Region Name=")
        RegionName = strtrim(strrep(line, "Region Name=", ""));
    elseif startsWith(line, "Lens Mode=")
        LensMode = strtrim(strrep(line, "Lens Mode=", ""));
    elseif startsWith(line, "Pass Energy=")
        PassEnergy = str2double(strtrim(strrep(line, "Pass Energy=", "")));
    elseif startsWith(line, "Number of Sweeps=")
        NumberOfSweeps = str2double(strtrim(strrep(line, "Number of Sweeps=", "")));
    elseif startsWith(line, "Excitation Energy=")
        ExcitationEnergy = str2double(strtrim(strrep(line, "Excitation Energy=", "").replace(",", ".")));
    elseif startsWith(line, "Energy Scale=")
        EnergyScale = strtrim(strrep(line, "Energy Scale=", ""));
    elseif startsWith(line, "Acquisition Mode=")
        AcquisitionMode = strtrim(strrep(line, "Acquisition Mode=", ""));
    elseif startsWith(line, "Energy Unit=")
        EnergyUnit = strtrim(strrep(line, "Energy Unit=", ""));
    elseif startsWith(line, "Step Time=")
        StepTime = str2double(strtrim(strrep(line, "Step Time=", "").replace(",", ".")));
    elseif startsWith(line, "Detector First X-Channel=")
        DetectorFirstXChannel = str2double(strtrim(strrep(line, "Detector First X-Channel=", "")));
    elseif startsWith(line, "Detector Last X-Channel=")
        DetectorLastXChannel = str2double(strtrim(strrep(line, "Detector Last X-Channel=", "")));
    elseif startsWith(line, "Detector First Y-Channel=")
        DetectorFirstYChannel = str2double(strtrim(strrep(line, "Detector First Y-Channel=", "")));
    elseif startsWith(line, "Detector Last Y-Channel=")
        DetectorLastYChannel = str2double(strtrim(strrep(line, "Detector Last Y-Channel=", "")));
    elseif startsWith(line, "Number of Slices=")
        NumberOfSlices = str2double(strtrim(strrep(line, "Number of Slices=", "")));
    % elseif startsWith(line, "File=")
    %     path = strtrim(strrep(line, "File=", "").replace("\\\\", "\\"));
    %     FilePath = path(1:strfind(path, '\', 1, 'last')) + fileName;
    elseif startsWith(line, "Sequence=")
        Sequence = strtrim(strrep(line, "Sequence=", "").replace("\\\\", "\\"));
    elseif startsWith(line, "Spectrum Name=")
        SpectrumName = strtrim(strrep(line, "Spectrum Name=", ""));
    elseif startsWith(line, "Instrument=")
        Instrument = strtrim(strrep(line, "Instrument=", ""));
    elseif startsWith(line, "Location=")
        Location = strtrim(strrep(line, "Location=", ""));
    elseif startsWith(line, "User=")
        User = strtrim(strrep(line, "User=", ""));
    elseif startsWith(line, "Sample=")
        Sample = strtrim(strrep(line, "Sample=", ""));
    elseif startsWith(line, "Comments=")
        Comments = strtrim(strrep(line, "Comments=", ""));
    elseif startsWith(line, "Date=")
        Date = strtrim(strrep(line, "Date=", ""));
    elseif startsWith(line, "Time=")
        Time = strtrim(strrep(line, "Time=", ""));
    elseif startsWith(line, "Time per Spectrum Channel=")
        TimePerSpectrumChannel = str2double(strtrim(strrep(line, "Time per Spectrum Channel=", "").replace(",", ".")));
    elseif startsWith(line, "DetectorMode=")
        DetectorMode = strtrim(strrep(line, "DetectorMode=", ""));
    elseif startsWith(line, "A=")
        ManipulatorAzimuth = str2double(strtrim(strrep(line, "A=", "").replace(",", ".")));
    elseif startsWith(line, "P=")
        ManipulatorPolar = str2double(strtrim(strrep(line, "P=", "").replace(",", ".")));
    elseif startsWith(line, "T=")
        ManipulatorTilt = str2double(strtrim(strrep(line, "T=", "").replace(",", ".")));
    elseif startsWith(line, "X=")
        ManipulatorX = str2double(strtrim(strrep(line, "X=", "").replace(",", ".")));
    elseif startsWith(line, "Y=")
        ManipulatorY = str2double(strtrim(strrep(line, "Y=", "").replace(",", ".")));
    elseif startsWith(line, "Z=")
        ManipulatorZ = str2double(strtrim(strrep(line, "Z=", "").replace(",", ".")));
    end
end

% Close the file
fclose(fileID);

% Load axis data
Energy = linspace(LowEnergy, ...
                        LowEnergy + (numEnergyPoints-1)*EnergyStep, ...
                        numEnergyPoints);

if startsWith(EnergyScale,'Kinetic')
   dataStr.AxisLabel{1} = 'Kinetic energy';
else
   dataStr.AxisLabel{1} = 'Binding energy'; 
end

dataStr.AxisUnits{1} = 'eV';

Angle = linspace(LowAnalyzerAngle, ...
                        LowAnalyzerAngle + (numAnalyzerPoints-1)*AnalyzerStep, ...
                        numAnalyzerPoints); 

dataStr.AxisLabel{2} = 'Analyzer angle';
dataStr.AxisUnits{2} = '$\degree$';

Tilt = linspace(LowDeflectorAngle, ...
                        LowDeflectorAngle + (numDeflectorPoints-1)*DeflectorStep, ...
                        numDeflectorPoints);

dataStr.AxisLabel{3} = 'Deflector angle';
dataStr.AxisUnits{3} = '$\degree$';                    
%% Load data

fid = fopen(spectrumFile_bin,'r');
if fid == -1
    error('File cannot be opened: %s', spectrumFile_bin);
end
binData = fread(fid,inf,'float32');
fclose(fid);

Data = zeros(numel(Energy),numel(Angle),numel(Tilt));
for i = 1:numel(Angle)
    idx = (i-1)*numel(Energy); 
    for j = 1:numel(Tilt)
       Data(:,i,j) = binData((j-1)*numel(Angle)*numel(Energy) + idx + (1:numel(Energy))); 
    end
end
    %% Assemble data file
    %meta info
    dataStr.FileName	= FileName;
    %dataStr.H5file 	= '';
    dataStr.TimeStamp   = [Date Time];
    dataStr.Type       = "Eb(kx,ky)";
    dataStr.index      = 1:size(Data, 3);
    dataStr.raw_data 	= Data;
    dataStr.raw_tht  	= Angle;
    dataStr.raw_eb   	= Energy';
    dataStr.hv         = ExcitationEnergy;
    dataStr.tltM     	= Tilt;
    dataStr.thtM       = ManipulatorTilt;
    dataStr.tltE        = ManipulatorPolar;
    dataStr.meta       = struct();
    dataStr.meta.info   = Comments;
    dataStr.meta.Epass  = PassEnergy;
    dataStr.meta.Mode   = LensMode;
    dataStr.meta.X      = ManipulatorX;
    dataStr.meta.Y      = ManipulatorY;
    dataStr.meta.Z      = ManipulatorZ;
    dataStr.meta.Azimuth= ManipulatorAzimuth;

    % Check if the directory exists before attempting to remove it
if isfolder(tempDir)
    rmdir(tempDir, 's'); % 's' option is used to remove all contents of the directory
else
    disp('Directory does not exist.');
end
end