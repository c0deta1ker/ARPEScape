function save_adress_data(dataStr, SaveFullName)
% save_adress_data(dataStr, SaveFullName, isoslice_args)
%   This function save the ARPES data that has been processed using the
%   PESTools package. The function saves the data as a MATLAB structure
%   file, which can be loaded back in at any point of the processing /
%   analysis stage. This function should only be used to save a single 
%   MATLAB structure (as a .mat file).
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   dataStr:            MATLAB data structure containing all the ADRESS data.
%   -   SaveFullName:       (if empty, it is prompted)  string or char of the full path + filename to be saved.
%   -   isoslice_args:   	(if empty, it is ignored)   1x4 cell of {scanIndex, isoType, isoWin, remap} for quickly plotting a snapshot of the data.
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
    filter = {'*.mat'};
    [save_filename, save_filepath] = uiputfile(filter, 'Save the data...', FileName);
    save_fullfile = char(string(save_filepath) + string(save_filename));
    % - If Cancel is pressed, then return nothing
    if isequal(save_filepath,0) || isequal(save_filename,0); return; end
else
    save_fullfile = char(string(SaveFullName) + ".mat");
end

%% 2 - Executing the saving of the data
wbar = waitbar(0.5,'Saving...'); 
save(char(save_fullfile), 'dataStr', '-v7.3');
disp('-> saved data : '); display(dataStr);
close(wbar);

%% Close wait-bar
waitbar(1,'Save complete!'); 
end