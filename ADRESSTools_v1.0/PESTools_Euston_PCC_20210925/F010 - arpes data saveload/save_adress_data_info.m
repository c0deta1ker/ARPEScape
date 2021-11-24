function save_adress_data_info(dataStr, SaveFullName)
% save_adress_data(dataStr, SaveFullName)
%   This function saves the ARPES data that has been processed using the
%   PESTools package. The function generates an output text-file of the 
%   ARPES data structure, so that the final values of analysis can easily
%   be seen/used.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   dataStr:            MATLAB data structure containing all the ADRESS data.
%   -   SaveFullName:       (if empty, it is prompted)  string or char of the full path + filename to be saved.
%
%   OUT: (none, only the file is saved)

%% Default parameters
if nargin < 2; SaveFullName = ''; end
if isempty(SaveFullName); SaveFullName = ''; end
% -- Extract the initial file name
if isfield(dataStr, 'FileName');    FileName = dataStr.FileName;
else;                               FileName = '';
end
%% 1 - User defined FileName and Path for the processed data
if isempty(SaveFullName)
    filter = {'*.txt'};
    [save_filename, save_filepath] = uiputfile(filter, 'Save the data...', FileName);
    save_fullfile = char(string(save_filepath) + string(save_filename));
    % - If Cancel is pressed, then return nothing
    if isequal(save_filepath,0) || isequal(save_filename,0); return; end
    
end
% - Verify the correct extension exists, otherwise add it
A = string(SaveFullName);
if A(1:end-4) ~= ".txt";    save_fullfile = char(string(SaveFullName));
else;                       save_fullfile = char(string(SaveFullName) + ".txt");
end

%% 2 - Saving a text file with all the info
diary_name = char(save_fullfile);
diary(diary_name);
display(dataStr)
display('hello there')
diary off;

%% Close wait-bar
waitbar(1,'Save complete!'); 
end