function dataStr = load_cassiopee_kz(PathName)
%% Default parameters
if nargin < 1; PathName = ''; end
% - Display text
disp('Loading SOLEIL data...')

% - Verifying inputs are characters
PathName = char(PathName);
fileExt = '*.txt';
dataStr = create_arpes_data();

% - If no filename if defined, load in an empty data structure
if isempty(PathName); return;
% - Verify that a .h5 file has been parsed and load the variables
else
    %% 1 - Loading and reading in the txt data file    
    % Define the regular expression pattern to match the numbers in the filenames
    pattern = '(\d+)_ROI1_\.txt';
    
    % Get a list of the files in the directory
    file_list = dir(fullfile(string(PathName), fileExt));
    
    % Extract the file numbers using regular expressions
    file_nums = zeros(size(file_list));
    for i = 1:numel(file_list)
        filename = file_list(i).name;
        matches = regexp(filename, pattern, 'tokens');
        if ~isempty(matches)
            file_nums(i) = str2double(matches{1});
        end
    end
    
    % Sort the file list based on the extracted numbers
    [sorted_nums, idx] = sort(file_nums);
    sortedFileNames = file_list(idx);

    % Loop through the sorted file names and extract the theta values
    data_cell = cell(length(sortedFileNames), 1);

    for i = 1:numel(sortedFileNames)
        file_path = fullfile(string(PathName), sortedFileNames(i).name);
        file_contents = fileread(file_path);
    
        data = readmatrix(file_path, 'Range', 'B46');
        data_cell{i} = data;
    
        deb_pattern = 'Energy Step=\s*(\S+)';
        deb_match = regexp(file_contents,deb_pattern,'tokens');
        deb = deb_match{1}{1};
    end
    
    dimension_1_scale_pattern = 'Dimension 1 scale=\s*([\d\. ]+)';
    dimension_1_scale_match = regexp(file_contents, dimension_1_scale_pattern, 'tokens');
    dimension_1_scale = str2num(dimension_1_scale_match{1}{1});
    
    dimension_1_name_pattern = 'Dimension 1 name=\s*(\S+)';
    dimension_1_name_match = regexp(file_contents, dimension_1_name_pattern, 'tokens');
    dimension_1_name = dimension_1_name_match{1}{1};
    
    dimension_2_name_pattern = 'Dimension 2 name=\s*(\S+)';
    dimension_2_name_match = regexp(file_contents, dimension_2_name_pattern, 'tokens');
    dimension_2_name = dimension_2_name_match{1}{1};
    
    temperature_pattern = 'inputA=\s*([\d\. ]+)';
    temperature_match = regexp(file_contents, temperature_pattern, 'tokens');
    Temp = str2double(temperature_match{1}{1});
    
    % Concatenate the matrices
    Data = cat(3, data_cell{:});
    
    % Read the last file for the analyzer angle
    fid = fopen(file_path, 'r');
    line_num = 0;
    while ~feof(fid)
        line_num = line_num + 1;
        line = fgetl(fid);
        if strncmp(line, 'Dimension 2 scale=', length('Dimension 2 scale='))
            % Found the line starting with "Dimension 2 scale="
            parts = strsplit(line, '=');
            C = textscan(parts{2},'%f');
            % extract the array
            Raw_tht = C{1};
            break;
        end
    end
    fclose(fid);
    
    
    %% Read info files

    % Define the info folder path
    Info_Path = fullfile(string(PathName),'info');

    % Define the regular expression pattern to match the numbers in the filenames
    pattern = '(\d+)_i\.txt';
    
    % Get a list of the files in the directory
    file_list = dir(fullfile(string(Info_Path), fileExt));
    
    % Extract the file numbers using regular expressions
    file_nums = zeros(size(file_list));
    for i = 1:numel(file_list)
        filename = file_list(i).name;
        matches = regexp(filename, pattern, 'tokens');
        if ~isempty(matches)
            file_nums(i) = str2double(matches{1});
        end
    end
    
    % Sort the file list based on the extracted numbers
    [sorted_nums, idx] = sort(file_nums);
    sortedFileNames = file_list(idx);
  
    % Loop through the sorted file names and extract the theta values
    theta_values = cell(numel(file_list), 1);
    hv_values = cell(numel(file_list), 1);
    
    for i = 1:numel(sortedFileNames)
        file_path = fullfile(Info_Path, sortedFileNames(i).name);
        file_contents = fileread(file_path);
    
        theta_matches = regexp(file_contents, 'theta \(deg\) :[\t\s]*([\d.-]+)', 'tokens');
        hv_matches = regexp(file_contents, 'hv \(eV\) :[ \t]*([\d.]+)', 'tokens');
    
        theta_values{i} = str2double(theta_matches{1}{1});
        hv_values{i} = str2double(hv_matches{1}{1});
    end
    
    
    %% Assemble data file
    %meta info
    dataStr.FileName	= sortedFileNames(1).name;
    %dataStr.H5file 	= '';
    dataStr.TimeStamp   = string(sortedFileNames(1).date);
    dataStr.Type       = "Eb(kx,kz)";
    dataStr.index      = 1:size(Data, 3);
    dataStr.meta       = struct();
    dataStr.raw_data 	= Data;
    dataStr.raw_tht  	= transpose(Raw_tht);
    dataStr.raw_eb   	= transpose(dimension_1_scale);
    dataStr.deb      	= str2double(deb);
    dataStr.hv         = transpose(cell2mat(hv_values));
    dataStr.dhv        = [];
    dataStr.tltM     	= theta_values{1};
    %dataStr.tltE     	= [];
    %dataStr.thtM       = [];
    dataStr.Temp   	= Temp;


end