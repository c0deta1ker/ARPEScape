function dataStr = load_sis_data(FileName, PathName)
% dataStr = load_sis_data(FileName, PathName)
%   This function loads in the HDF5 data files of ARPES data from 
%   the SIS beamline at the SLS. The output is a MATLAB
%   data-structure that yields all the data and information. This function 
%   should be used to load in a single SIS data-file that is unprocessed
%   or has been previously processed using PESTools.
%   REQ. FUNCTIONS:
%   -   [data, axes, note]  = ReaderHDF5(fname);
%   -   [DataXC,OffsE]      = SumScanXC(Energy,Data,maxLagE[,WinE]);
%
%   IN:
%   -   FileName:           char of the input .h5 or .mat file-name.
%   -   PathName:           char of the input .h5 or .mat directory path.
%
%   OUT:
%   dataStr - MATLAB data structure containing all fields below;
%	 .(FileName):        string of the current filename of the data file.
%	 .(H5file):          string of the raw .H5 filename of the data file.
%	 .(Type):            string "Eb(k)", "Eb(k,i)", "Eb(kx,ky)" or "Eb(kx,kz)".
%	 .(raw_data):        2Data [nEnergy x nAngle] or 3Data [nEnergy x nAngle xnScan] data array.
%	 .(raw_tht):         [1×nAngle] row vector of angles.
% 	 .(raw_eb):          [nEnergyx1] column vector of energy.
%	 .(hv):              scalar or [1×nScan] row vector (if Scan parameter for 3Data data).
%	 .(tltM):            scalar or [1×nScan] row vector (if Scan parameter for 3Data data).
%	 .(thtM):            scalar of the manipular Theta angle.
%	 .(Temp):            scalar of the sample temperature.
%	 .(meta):            structure with all measurement information.

%% Default parameters
if nargin < 2; PathName = '';  end
if isempty(PathName);   PathName = '';  end
disp('Loading in SIS data...')
% - Verifying input are characters
FileName    = char(FileName);
PathName    = char(PathName);

%% 1 - Reading in the HDataF5 data
% - Verify that a .h5 file has been parsed and load the variables
if string(FileName(end-2:end)) == ".h5"
    % -- Defining the MATLAB data-structure
    dataStr             = struct();
    dataStr.FileName    = FileName(1:end-3);
    dataStr.H5path      = char(PathName);
    dataStr.H5file      = char(FileName);
    H5full              = char(string(PathName) + string(FileName));
    
    %% -- Extracting all of the data variables
    datav_new           = struct();
    datav_new.value     = h5read(H5full,'/Electron Analyzer/Image Data');
    if ndims(datav_new.value) == 2
        size3       = size(datav_new.value);
        xnum        = size3(1);ynum=size3(2);
        xaxis       = h5readatt(H5full,'/Electron Analyzer/Image Data','Axis1.Scale');
        datav_new.y = xaxis(1):xaxis(2):xaxis(1)+(xnum-1)*xaxis(2);
        yaxis       = h5readatt(H5full,'/Electron Analyzer/Image Data','Axis0.Scale');
        datav_new.z = yaxis(1):yaxis(2):yaxis(1)+(ynum-1)*yaxis(2);
    elseif ndims(datav_new.value) == 3
        size3       = size(datav_new.value);
        xnum        = size3(1);ynum=size3(3);znum=size3(2);
        xaxis       = h5readatt(H5full,'/Electron Analyzer/Image Data','Axis2.Scale');
        datav_new.x = xaxis(1):xaxis(2):xaxis(1)+(xnum-1)*xaxis(2);
        yaxis       = h5readatt(H5full,'/Electron Analyzer/Image Data','Axis0.Scale');
        datav_new.z = yaxis(1):yaxis(2):yaxis(1)+(ynum-1)*yaxis(2);
        zaxis       = h5readatt(H5full,'/Electron Analyzer/Image Data','Axis1.Scale');
        datav_new.y = zaxis(1):zaxis(2):zaxis(1)+(znum-1)*zaxis(2);
    end

    %% -- Reading in all of the meta data
    for meta = {h5info(H5full ,'/Electron Analyzer/').Datasets.Attributes.Name}
        meta= meta(1);
        if strcmp(meta{1},'Acquisition Mode')
            datav_new.info.Aquisition_Mode= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Acquisition Mode');
        end
        if strcmp(meta{1}, 'Intensity Units')
            datav_new.info.Intensity_Units= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Intensity Units');
        end
        if strcmp(meta{1}, 'Axis0.ScaleType')
            datav_new.info.Axis0.ScaleType= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Axis0.ScaleType');
        end
        if strcmp(meta{1}, 'Axis0.Scale')
            datav_new.info.Axis0.Scale= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Axis0.Scale');
        end
        if strcmp(meta{1}, 'Axis0.Description')
            datav_new.info.Axis0.Description= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Axis0.Description');
        end
        if strcmp(meta{1},'Axis0.Units')
            datav_new.info.Axis0.Units= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Axis0.Units');
        end
        if strcmp(meta{1},'Axis0.Mode')
            datav_new.info.Axis0.Mode= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Axis0.Mode');
        end
        if strcmp(meta{1}, 'Axis1.ScaleType')
            datav_new.info.Axis1.ScaleType= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Axis1.ScaleType');
        end
        if strcmp(meta{1}, 'Axis1.Scale')
            datav_new.info.Axis1.Scale= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Axis1.Scale');
        end
        if strcmp(meta{1}, 'Axis1.Description')
            datav_new.info.Axis1.Description= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Axis1.Description');
        end
        if strcmp(meta{1},'Axis1.Units')
            datav_new.info.Axis1.Units= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Axis1.Units');
        end
        if strcmp(meta{1},'Axis1.Mode')
            datav_new.info.Axis1.Mode= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Axis1.Mode');
        end
        if strcmp(meta{1},'Specified Number of Sweeps')
            datav_new.info.Sweeps= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Specified Number of Sweeps');
        end
        if strcmp(meta{1},'Dither Window (%)')
            datav_new.info.Dither.Window= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Dither Window (%)');
        end
        if strcmp(meta{1},'Dither Points')
            datav_new.info.Dither.Points= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Dither Points');
        end
        if strcmp(meta{1},'Lens Mode')
            datav_new.info.Lens_Mode= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Lens Mode');
        end
        if strcmp(meta{1},'Pass Energy (eV)')
            datav_new.info.EPass= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Pass Energy (eV)');
        end
        if strcmp(meta{1},'Energy Scale')
            datav_new.info.Energy_Scale= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Energy Scale');
        end
        if strcmp(meta{1},'Work Function (eV)')
            datav_new.info.Work_Function= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Work Function (eV)');
        end
        if strcmp(meta{1},'Dwell Time (ms)')
            datav_new.info.Dwell_Time= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Dwell Time (ms)');
        end
        if strcmp(meta{1},'Detector First X-Channel')
            datav_new.info.Detector_First_X_Channel= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Detector First X-Channel');
        end
        if strcmp(meta{1},'Detector Last X-Channel')
            datav_new.info.Detector_Last_X_Channel= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Detector Last X-Channel');
        end
        if strcmp(meta{1},'Detector First Y-Channel')
            datav_new.info.Detector_First_Y_Channel= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Detector First Y-Channel');
        end
        if strcmp(meta{1},'Detector Last Y-Channel')
            datav_new.info.Detector_Last_Y_Channel= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Detector Last Y-Channel');
        end
        if strcmp(meta{1},'Detector Slices')
            datav_new.info.Detector_Slices= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Detector Slices');
        end
        if strcmp(meta{1},'Detector ADC Mode')
            datav_new.info.Detector_ADC_Mode= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Detector ADC Mode');
        end
        if strcmp(meta{1},'Date Created')
            datav_new.info.Creation_Date= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Date Created');
        end
        if strcmp(meta{1},'Time Created')
            datav_new.info.Creation_Time= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Time Created');
        end
        if strcmp(meta{1},'Excitation Energy (eV)')
            datav_new.info.Excitation_Energy= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Excitation Energy (eV)');
        end
        if strcmp(meta{1},'Excitation Energy Read Device')
            datav_new.info.Excitation_Energy_Readback= ...
              h5readatt(H5full,'/Electron Analyzer/Image Data',...
               'Excitation Energy Read Device');
        end
    end
    for meta = {h5info(H5full ,'/').Attributes.Name}
        meta = meta(1);
        if strcmp(meta{1},'Comments')
            datav_new.info.Comments= ...
              h5readatt(H5full,'/',...
              'Comments');
        end
    end
    %% -- Reading in all of the instrumental data
    try datav_new.info.ExitSlit=h5read(H5full,...
            '/Other Instruments/Exit Slit'); catch ;end 
    try datav_new.info.FE_Vert_Width=h5read(H5full,...
            '/Other Instruments/FE Vert. Width'); catch ;end 
    try datav_new.info.FE_Horiz_Width=h5read(H5full,...
            '/Other Instruments/FE Horiz. Width'); catch ;end 
    try datav_new.info.Phi=h5read(H5full,...
            '/Other Instruments/Phi'); catch ;end 
    try datav_new.info.Pressure_AC=h5read(H5full,...
            '/Other Instruments/Pressure AC (ACMI)'); catch ;end 
    try datav_new.info.Pressure_TC=h5read(H5full,...
            '/Other Instruments/Pressure TC (TCMP)'); catch ;end 
    try datav_new.info.Ring_Current=h5read(H5full,...
            '/Other Instruments/Storage Ring Current'); catch ;end 
    try datav_new.info.Temperature_Cryo=h5read(H5full,...
            '/Other Instruments/Temperature A (Cryostat)'); catch ;end 
    try datav_new.info.Temperature_Sample1=h5read(H5full,...
            '/Other Instruments/Temperature B (Sample 1)'); catch ;end 
    try datav_new.info.Temperature_Head=h5read(H5full,...
            '/Other Instruments/Temperature C (Head Mech)'); catch ;end 
    try datav_new.info.Temperature_Sample2=h5read(H5full,...
            '/Other Instruments/Temperature D (Sample 2)'); catch ;end 
    try datav_new.info.Temperature_Boot1=h5read(H5full,...
            '/Other Instruments/Temperature D2 (Boot 1)'); catch ;end 
    try datav_new.info.Temperature_Shield=h5read(H5full,...
            '/Other Instruments/Temperature D3 (Shield)'); catch ;end 
    try datav_new.info.Temperature_Boot2=h5read(H5full,...
            '/Other Instruments/Temperature D4 (Boot 2)'); catch ;end 
    try datav_new.info.Temperature_Cryopump=h5read(H5full,...
            '/Other Instruments/Temperature D5 (Cryopump)'); catch ;end 
    try datav_new.info.Theta=h5read(H5full,...
            '/Other Instruments/Theta'); catch ;end
    try datav_new.info.Tilt=h5read(H5full,...
            '/Other Instruments/Tilt'); catch ;end 
    try datav_new.info.X=h5read(H5full,...
            '/Other Instruments/X'); catch ;end 
    try datav_new.info.Y=h5read(H5full,...
            '/Other Instruments/Y'); catch ;end 
    try datav_new.info.Z=h5read(H5full,...
            '/Other Instruments/Z'); catch ;end 
    try datav_new.info.hv=h5read(H5full,...
            '/Other Instruments/hv'); catch ;end 
    
    %% -- Creating a consistent data structure to be used
    % Checking to see whether the Scan parameter is identical
    if length(datav_new.info.Tilt) > 1 && all(datav_new.info.Tilt == datav_new.info.Tilt(1))
        index              	= 1:length(datav_new.info.Tilt);
        datav_new.info.Tilt = datav_new.info.Tilt(1);
    elseif length(datav_new.info.hv) > 1 && all(datav_new.info.hv == datav_new.info.hv(1))
        index               = 1:length(datav_new.info.hv);
        datav_new.info.hv   = datav_new.info.hv(1);
    elseif length(datav_new.info.Tilt) > 1
        index = 1:length(datav_new.info.Tilt);
    elseif length(datav_new.info.hv) > 1
        index = 1:length(datav_new.info.hv);
    else
        index = 1;
    end
    % Extracting the type of scan that is performed
    if length(datav_new.info.Tilt) > 1;     dataStr.Type = "Eb(kx,ky)";
    elseif length(datav_new.info.hv) > 1;   dataStr.Type = "Eb(kx,kz)";
    elseif length(index) > 1;               dataStr.Type = "Eb(k,i)";
    else;                                   dataStr.Type = "Eb(k)";
    end
    % Squeezingthe data if the scan parameter is 1 unit
    if size(datav_new.value, 1) == 1; datav_new.value = squeeze(datav_new.value); end
    
    % Assigning all data variables
    dataStr.index       = index;
    dataStr.raw_data    = datav_new.value;
    dataStr.raw_tht     = datav_new.y;
    dataStr.raw_eb      = datav_new.z;
    dataStr.hv          = datav_new.info.hv;
    dataStr.tltM        = datav_new.info.Tilt;
    dataStr.thtM        = datav_new.info.Theta;
    dataStr.Temp        = datav_new.info.Temperature_Sample1;
    dataStr.meta        = datav_new.info;
    
    % Validity checks on the data variables
    if dataStr.index == 1;  dataStr.raw_data = permute(dataStr.raw_data, [2 1]);
    else;                   dataStr.raw_data = permute(dataStr.raw_data, [3 2 1]);
    end
    if size(dataStr.raw_eb, 2) > 1;     dataStr.raw_eb = dataStr.raw_eb'; end
    if size(dataStr.raw_tht, 1) > 1;    dataStr.raw_tht = dataStr.raw_tht'; end
    if size(dataStr.tltM, 1) > 1;       dataStr.tltM = dataStr.tltM'; end
    if size(dataStr.hv, 1) > 1;         dataStr.hv = dataStr.hv'; end
    
%% 2 - Reading in the .mat data that has been processed
elseif string(FileName(end-3:end)) == ".mat"
    arpes_data  = load(char(string(PathName) + string(FileName)));
    dataStr     = arpes_data.dataStr;
    if isstruct(dataStr) && isfield(dataStr, 'FileName')
        dataStr.FileName = FileName(1:end-4);
    end
end

end