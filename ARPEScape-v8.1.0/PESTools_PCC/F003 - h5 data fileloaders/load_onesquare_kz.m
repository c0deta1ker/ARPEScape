function dataStr = load_onesquare_kz(PathName)
%% Default parameters
if nargin < 1; PathName = ''; end
% - Display text
disp('Loading One-Square data...')

% - Verifying inputs are characters
PathName = char(PathName);
fileExt = '*.ibw';
dataStr = create_arpes_data();

% - If no filename if defined, load in an empty data structure
if isempty(PathName); return;
else
    %% Loading and reading the data files    
   
    % Get a list of the files in the directory
    file_list = dir(fullfile(string(PathName), fileExt));
    
    file_nums = zeros(size(file_list));
  
    % Sort the file list based on the extracted numbers
    [sorted_nums, idx] = sort(file_nums);
    sortedFileNames = file_list(idx);

    % Read first scan information
    SliceData = IBWread(sortedFileNames(1).name, string(PathName));
    thtM = SliceData.thtM;
    tltE = SliceData.tltE;
    tltM = SliceData.tltM;
    raw_tht = SliceData.raw_tht;
    raw_eb = SliceData.raw_eb;

    % Loop through the file list
    data_cell = cell(numel(file_list), 1);
    hv_values = cell(numel(file_list), 1);

    for i = 1:numel(sortedFileNames)
        SliceData = IBWread(sortedFileNames(i).name, string(PathName));
        
        data_cell{i} = SliceData.raw_data;
        hv_values{i} = SliceData.hv;
    end
    
    % Concatenate the matrices
    Data = cat(3, data_cell{:});
    
    %% Assemble data file
    %meta info
    dataStr.FileName	= sortedFileNames(1).name;
    %dataStr.H5file 	= '';
    dataStr.TimeStamp   = string(sortedFileNames(1).date);
    dataStr.Type        = "Eb(kx,kz)";
    dataStr.index       = 1:size(Data, 3);
    % - Assigning the meta data
    %dataStr.meta.deb    = str2double(deb);
    %dataStr.meta.dhv    = [];
    % - Assinging experimental parameters
    dataStr.hv          = transpose(cell2mat(hv_values));
    dataStr.tltM     	= tltM;
    dataStr.tltE     	= tltE;
    dataStr.thtM        = thtM;
    %dataStr.Temp   	    = Temp;
    % - Assinging the main experimental variables
    dataStr.raw_data 	= Data;
    dataStr.raw_tht  	= raw_tht;
    dataStr.raw_eb   	= raw_eb;


end