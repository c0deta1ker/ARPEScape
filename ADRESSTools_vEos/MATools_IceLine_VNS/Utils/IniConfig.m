classdef IniConfig < handle
    %IniConfig - The class for working with configurations of settings and INI-files. 
    % This class allows you to create configurations of settings, and to manage them. 
    % The structure of the storage settings is similar to the structure of 
    % the storage the settings in the INI-file format. 
    % The class allows you to import settings from the INI-file and to export 
    % the settings in INI-file. 
    % Can be used for reading/writing data in the INI-file and managing 
    % settings of application.
    %
    %
    % Config Syntax:
    %
    %   ; some comments
    %
    %   [Section1] ; allowed the comment to section
    %   ; comment on the section
    %   key1 = value_1 ; allowed a comment to an individual parameter
    %   key2 = value_2
    %   
    %   [Section2]
    %   key1 = value_1.1, value_1.2     ; array data
    %   key2 = value_2
    %   ...
    %
    % Note:
    %   * There may be spaces in the names of sections and keys
    %   * Keys should not be repeated in the section (will read the last)
    %   * Sections should not be repeated in the config (will read the last)
    % ---------------------------------------------------------------------
    %
    % Support data types:
    %   * numeric scalars and vectors
    %   * logical scalars and vectors
    %   * strings
    %
    % ---------------------------------------------------------------------
    %
    % Using:
    %   ini = IniConfig();  % Constructor (an instance of class)
    %
    % Public Properties:
    %   Enter this command to get the properties:
    %   >> properties IniConfig
    %
    % Public Methods:
    %   Enter this command to get the methods:
    %   >> methods IniConfig
    %
    %   Enter this command to get more info of method:
    %   >> help IniConfig/methodname
    %
    % ---------------------------------------------------------------------
    %
    % Example:
    %   ini = IniConfig();
    %   ini.ReadFile('example.ini')
    %   sections = ini.GetSections()
    %   [keys, count_keys] = ini.GetKeys(sections{1})
    %   values = ini.GetValues(sections{1}, keys, cell(count_keys, 1))
    %   new_values = cell(count_keys, 1);
    %   new_values(:) = {rand()};
    %   ini.SetValues(sections{1}, keys, new_values)
    %   ini.WriteFile('example1.ini')
    %
    %
    % Example:
    %   ini = IniConfig();
    %   ini.AddSection('Some Section 1')
    %   ini.AddKey('Some Section 1', 'some_key1', 'hello!')
    %   ini.AddKey('Some Section 1', 'some_key2', [10, 20])
    %   ini.AddKey('Some Section 1', 'some_key3', true)
    %   ini.AddSection('Some Section 2')
    %   ini.AddKey('Some Section 2', 'some_key1')
    %   ini.WriteFile('example2.ini')
    %
    %
    % Example:
    %   ini = IniConfig();
    %   ini.AddSection('Some Section 1')
    %   ini.AddKey('Some Section 1', 'some_key1', 'hello!')
    %   ini.AddKey('Some Section 1', 'some_key2', [10, 20])
    %   ini.AddKey('Some Section 1', 'some_key3', [false, true])
    %   ini.WriteFile('example31.ini')
    %   ini.RemoveKey('Some Section 1', 'some_key1')
    %   ini.RenameKey('Some Section 1', 'some_key2', 'renamed_some_key2')
    %   ini.RenameSection('Some Section 1', 'Renamed Section 1')
    %   ini.WriteFile('example32.ini')
    %
    %
    % ---------------------------------------------------------------------
    % Note:
    %   This implementation does not use cycles, and works fairly quickly,
    %   even on large files. 
    %   Parsing the files by using regular expressions
    %   All comments are saved when set values
    %
    %   This class was tested on the win32 platform only but it should
    %   also work on Unix/Linux platforms.
    %
    %
    % ---------------------------------------------------------------------
    % Author:         Iroln <iroln85@gmail.com>
    % Version:        0.6.2
    % First release:  25.07.09
    % Last revision:  11.08.09
    % Copyright:      (c) 2009 Evgeny Pr aka Iroln. All rights reserved.
    %
    % Bug reports, questions, etc. can be sent to the e-mail given above.
    % ---------------------------------------------------------------------
    
    
    properties (GetAccess = 'public', SetAccess = 'private')
        comment_style = ';'
        count_sections = 0
        full_count_keys = 0
    end
    
    properties (GetAccess = 'private', SetAccess = 'private')
        file_name = ''
        config_data_array = {}
        indexes_of_sections = []
        indexes_of_empty_str = []
        
        count_strings = 0
        count_empty_strings = 0
        
        is_created_configuration = false
    end
    
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Public Methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %------------------------------------------------------------------
        function obj = IniConfig()
            %IniConfig - Constructor
            % To Create empty default configuration.
            % Using:
            %   ini_obj = IniConfig()
            %
            
            obj.CreateIni();
        end
        
        %------------------------------------------------------------------
        function status = ReadStream(obj, data_stream, comment_style)
            obj.comment_style = ValidateCommentStyle(comment_style);
            file_data = data_stream;
            obj.count_strings = size(file_data, 1);
                
            % Parse data
            obj.config_data_array = IniFileParser(file_data, obj.comment_style);    
            obj.UpdateSectionsInfo()
            obj.UpdateEmptyStringsInfo()
            obj.UpdateCountKeysInfo()    

            obj.is_created_configuration = true;            
        end
        
        %------------------------------------------------------------------
        function status = ReadFile(obj, filename, comment_style)
            %ReadFile - to read in the object the config data from a INI file
            % Using:
            %   status = ReadFile(filename, comment_style)
            %       status - 1 (true): Success
            %       status - 0 (false): Failed
            %
            
            if (nargin < 2)
                error('Not enough input arguments.')
            elseif (nargin == 3)
                obj.comment_style = ValidateCommentStyle(comment_style);
            end
            
            if ~ischar(filename)
                error('Requires string input.')
            else
                % Get data from file
                fid = fopen(filename, 'r');
                
                if (fid ~= -1)
                    file_data = textscan(fid, ...
                        '%s', ...
                        'delimiter', '\n', ...
                        'endOfLine', '\r\n');
                    
                    fclose(fid);
                    status = true;
                else
                    status = false;
                    return;
                end
                
                obj.file_name = filename;
                file_data = file_data{1};
                obj.count_strings = size(file_data, 1);
                
                % Parse data
                obj.config_data_array = IniFileParser(file_data, obj.comment_style);
                
                obj.UpdateSectionsInfo()
                obj.UpdateEmptyStringsInfo()
                obj.UpdateCountKeysInfo()
                
                obj.is_created_configuration = true;
            end
        end
        
        %------------------------------------------------------------------
        function status = IsSection(obj, section_name)
            %IsSection - determine whether there is a section
            % Using:
            %   status = IsSection(section_name)
            %       status - 1 (true): section is found
            %       status - 0 (false): section is not found
            %
            
            if (obj.is_created_configuration)
                if (nargin < 2)
                    error('Not enough input arguments.')
                end
                
                section_name = obj.ValidateSectionName(section_name);
                sect_num = obj.NameToNumSection(section_name);
                
                if ~isempty(sect_num)
                    status = ~obj.IsUniqueSectionName(section_name);
                else
                    status = false;
                end
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function status = IsKey(obj, section_name, key_name)
            %IsSection - determine whether there is a key in a given section
            % Using:
            %   status = IsKey(section_name, key_name)
            %       status - 1 (true): key is found
            %       status - 0 (false): key is not found
            %
            
            if (obj.is_created_configuration)
                if (nargin < 3)
                    error('Not enough input arguments.')
                end
                
                section_name = obj.ValidateSectionName(section_name);
                key_name = obj.ValidateKeyName(key_name);
                sect_num = obj.NameToNumSection(section_name);
                
                if (~isempty(sect_num) && ~isempty(key_name))
                    status = ~obj.IsUniqueKeyName(section_name, key_name);
                else
                    status = false;
                end
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function [names_sect, count_sect] = GetSections(obj)
            %GetSections - get names of all sections
            % Using:
            %   [names_sect, count_sect] = GetSections()
            %
            % Output:
            %   names_sect - cell array with the names of sections
            %   count_sect - number of sections in configuration
            %
            
            if (obj.is_created_configuration)
                names_sect = obj.config_data_array(obj.indexes_of_sections, 1);
                count_sect = obj.count_sections;
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function [keys_names, count_keys] = GetKeys(obj, section_name)
            %GetKeys - get names of all keys from given section
            % Using:
            %   [keys_names, count_keys] = GetKeys(section_name)
            %
            % Output:
            %   keys_names - cell array with the names of keys
            %   count_keys - number of keys in given section
            %
            
            if (obj.is_created_configuration)
                if (nargin < 2)
                    error('Not enough input arguments.')
                end
                
                section_name = obj.ValidateSectionName(section_name);
                
                [keys_indexes, count_keys] = obj.GetKeysIndexes(section_name);
                
                if isempty(keys_indexes)
                    keys_names = {};
                    count_keys = 0;
                    return;
                end
                
                keys_names = obj.config_data_array(keys_indexes, 1);
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function [values, count_read] = GetValues(obj, section_name, keys_names, default_values, param_name, param_val)
            %GetValues - get values of keys from given section
            % Using:
            %   [values, count_read] = GetValues(section_name, ...
            %                                    keys_names, ...
            %                                    default_values, ...
            %                                    param_name, ...
            %                                    param_val)
            % Output:
            %   values - cell array with the values
            %   count_read - number of successfully get values
            %
            % Input:
            %   section_name - name of given section
            %   keys_names - names of given keys
            %   default_values - default values
            %
            % Optional Input Parameters:
            %   param_name - must be 'AlwaysReturnCell'
            %   param_val - must be 0 or 1
            %
            
            if (obj.is_created_configuration)
                if (nargin < 4)
                    error('Not enough input arguments.')
                elseif (nargin == 4)
                    is_get_values_always_in_cell = true;
                    
                elseif (nargin == 5)
                    error('Not enough input arguments.')
                    
                elseif (nargin == 6)
                    if ischar(param_name)
                        if strcmpi('AlwaysReturnCell', param_name)
                            if (param_val == 0 || param_val == 1)
                                is_get_values_always_in_cell = logical(param_val);
                            else
                                error('Value must be 0 or 1.')
                            end
                        else
                            error('Invalid param name: "%s". Must be "AlwaysReturnCell"', param_name)
                        end
                    else
                        error('The method GetValue expects a string "AlwaysReturnCell".')
                    end
                end
                
                section_name = obj.ValidateSectionName(section_name);
                
                [keys_names, default_values] = ...
                    ValidateInputsForGetAndSetValues(section_name, ...
                    keys_names, default_values);
                
                count_keys = numel(keys_names);
                values = cell(count_keys, 1);
                count_read = 0;
                
                keys_indexes = obj.GetKeysIndexes(section_name);
                
                if isempty(keys_indexes)
                    values = default_values;
                    
                    if (length(values) == 1 && ~is_get_values_always_in_cell)
                        values = values{1};
                    end
                    
                    return;
                end
                
                keys_names = strtrim(keys_names);
                keys = obj.config_data_array(keys_indexes, 1);
                vals = obj.config_data_array(keys_indexes, 2);
                
                % Get the names of keys, which are not in the INI file 
                % in a given section and specifying this key Defaults
                [not_found_keys, not_found_indexes] = setdiff(keys_names, keys);
                values(not_found_indexes) = default_values(not_found_indexes);
                
                % Getting the names of keys, which is in the INI file in a 
                % given section
                [found_keys, found_indexes, found_data_indexes] = intersect(keys_names, keys);
                values(found_indexes) = vals(found_data_indexes);
                count_read = length(found_indexes);
                
                %DEBUG: found_indexes = sort(found_indexes);
                values(found_indexes) = ParserValues(values(found_indexes));
                
                if (length(values) == 1 && ~is_get_values_always_in_cell)
                    values = values{1};
                end
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function count_write = SetValues(obj, section_name, keys_names, values)
            %SetValues - set values for given keys from given section
            % Using:
            %   count_write = SetValues(section_name, keys_names, values)
            %
            % Output:
            %   count_write - number of successfully set values
            %
            % Input:
            %   section_name - name of given section (must be string)
            %   keys_names - names of given keys (must be cell array of strings or string)
            %   values - values of keys (must be cell array or one value)
            %
            % Note:
            %   Number of elements in the second and third input argument must be equal
            %
            
            if (obj.is_created_configuration)
                if (nargin < 4)
                    error('Not enough input arguments.')
                end
                
                section_name = obj.ValidateSectionName(section_name);
                
                [keys_names, values] = ...
                    ValidateInputsForGetAndSetValues(section_name, ...
                    keys_names, values);
                
                count_keys = numel(keys_names);
                count_write = 0;
                keys_indexes = obj.GetKeysIndexes(section_name);
                
                if isempty(keys_indexes)
                    return;
                end
                
                [valid_keys, count_valid_keys] = obj.GetKeys(section_name);
                
                if (count_keys > count_valid_keys)
                    values = values(1:count_valid_keys);
                end
                
                keys_names = strtrim(keys_names);
                keys = obj.config_data_array(keys_indexes, 1);
                
                % If the key is not unique within the section, it will be 
                % processed by the key, which is located the last
                [keys_names, m] = unique(keys_names);
                values = values(m);
                
                [not_found_keys, not_found_indexes] = setdiff(keys_names, keys);
                keys_names(not_found_indexes) = [];
                values(not_found_indexes) = [];
                
                [found_keys, found_indexes, found_data_indexes] = intersect(keys_names, keys);
                values(sort(found_indexes)) = values(found_indexes);
                
                % Get string value indices
                str_indexes = cellfun(@(x) ischar(x), values);
                
                % Get numeric value indices
                num_indexes = cellfun(@(x) isnumeric(x), values);
                
                % Get logical value indices
                log_indexes = cellfun(@(x) islogical(x), values);
                
                values(num_indexes | log_indexes) = cellfun(@(x) mat2str(x), ...
                    values(num_indexes | log_indexes), 'UniformOutput', 0);
                
                values(~str_indexes & ~num_indexes & ~log_indexes) = [];
                values(str_indexes) = strtrim(values(str_indexes));
                
                values(num_indexes | log_indexes) = ...
                    CorrectNumArrayStrings(values(num_indexes | log_indexes));
                
                if isempty(values)
                    values = {''};
                end
                
                obj.config_data_array(keys_indexes(found_data_indexes), 2) = values;
                count_write = length(values);
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function CreateIni(obj)
            %CreateIni - create new empty INI configuration
            % Using:
            %   CreateIni()
            %
            
            obj.config_data_array = cell(2,3);
            obj.config_data_array(:,:) = {''};
            
            obj.UpdateCountStrings();
            obj.UpdateSectionsInfo()
            obj.UpdateEmptyStringsInfo()
            
            obj.file_name = '';
            obj.is_created_configuration = true;
        end
        
        %------------------------------------------------------------------
        function status = AddSection(obj, section_name)
            %AddSection - add section to end configuration
            % Using:
            %   status = AddSection(section_name)
            %       status - 1 (true): Success
            %       status - 0 (false): Failed
            %
            
            if (obj.is_created_configuration)
                if (nargin < 2)
                    error('Not enough input arguments.')
                end
                
                section_name = obj.ValidateSectionName(section_name);
                status = obj.InsertSection(obj.count_sections+1, section_name);
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function status = InsertSection(obj, pos, section_name)
            %InsertSection - insert section to given position
            % Using:
            %   status = InsertSection(pos, section_name)
            %       status - 1 (true): Success
            %       status - 0 (false): Failed
            %
            
            if (obj.is_created_configuration)
                if (nargin < 3)
                    error('Not enough input arguments.')
                else
                    if (~isnumeric(pos) && ~isscalar(pos))
                        error('Position must be a scalar positive integer')
                    end
                    
                    if (pos < 1)
                        error('Position must be a positive integer')
                    elseif (pos > obj.count_sections+1)
                        pos = obj.count_sections+1;
                    end
                end
                
                section_name = obj.ValidateSectionName(section_name);
                
                if ~isempty(section_name)
                    is_unique_sect = obj.IsUniqueSectionName(section_name);
                    if ~is_unique_sect
                        status = false;
                        return;
                    end

                    if (pos <= obj.count_sections && obj.count_sections > 0)
                        sect_ind = obj.indexes_of_sections(pos);
                        
                    elseif (pos == 1 && obj.count_sections == 0)
                        sect_ind = 1;
                        obj.config_data_array = {};
                        
                    elseif (pos == obj.count_sections+1)
                        %FIXME:
                        sect_ind = obj.count_strings+1;
                    end
                    
                    new_data = cell(2,3);
                    new_data(1,:) = {section_name, '', ''};
                    new_data(2,:) = {''};
                    
                    obj.config_data_array = insertcell(obj.config_data_array, ...
                        sect_ind, new_data);
                    
                    obj.UpdateCountStrings();
                    obj.UpdateSectionsInfo();
                    obj.UpdateEmptyStringsInfo();
                    
                    status = true;
                else
                    status = false;
                end
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function status = RemoveSection(obj, section_name)
            %RemoveSection - remove given section
            % Using:
            %   status = RemoveSection(section_name)
            %       status - 1 (true): Success
            %       status - 0 (false): Failed
            %
            
            if (obj.is_created_configuration)
                if (nargin < 2)
                    error('Not enough input arguments.')
                end
                
                section_name = obj.ValidateSectionName(section_name);
                sect_num = obj.NameToNumSection(section_name);
                
                if ~isempty(sect_num)
                    
                    if (sect_num < obj.count_sections)
                        first_ind = obj.indexes_of_sections(sect_num);
                        last_ind = obj.indexes_of_sections(sect_num+1)-1;
                        
                    elseif (sect_num == obj.count_sections)
                        first_ind = obj.indexes_of_sections(sect_num);
                        last_ind = obj.count_strings;
                    end
                    
                    obj.config_data_array(first_ind:last_ind,:) = [];
                    
                    obj.UpdateCountStrings();
                    obj.UpdateSectionsInfo();
                    obj.UpdateEmptyStringsInfo();
                    obj.UpdateCountKeysInfo();
                    
                    status = true;
                else
                    status = false;
                end
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function status = RenameSection(obj, section_name, new_section_name)
            %RenameSection - rename given section
            % Using:
            %   status = RenameSection(section_name, new_section_name)
            %       status - 1 (true): Success
            %       status - 0 (false): Failed
            %
            
            if (obj.is_created_configuration)
                if (nargin < 3)
                    error('Not enough input arguments.')
                end
                
                section_name = obj.ValidateSectionName(section_name);
                new_section_name = obj.ValidateSectionName(new_section_name);
                sect_num = obj.NameToNumSection(section_name);
                
                if (~isempty(new_section_name) && ~isempty(sect_num))
                    sect_ind = obj.indexes_of_sections(sect_num);
                    
                    obj.config_data_array(sect_ind, 1) = {new_section_name};
                    status = true;
                else
                    status = false;
                end
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function [status, write_value] = AddKey(obj, section_name, key_name, value)
            %AddKey - add key in a end given section
            % Using:
            %   [status, write_value] = AddKey(section_name, key_name, value)
            %
            % Output:
            %   status - 1 (true): Success, status - 0 (false): Failed
            %   write_value - 1 (true): Success, status - 0 (false): Failed
            %
            
            if (obj.is_created_configuration)
                if (nargin < 3)
                    error('Not enough input arguments.')
                elseif (nargin == 3)
                    value = '';
                end
                
                section_name = obj.ValidateSectionName(section_name);
                [keys_indexes, count_keys] = obj.GetKeysIndexes(section_name);
                
                [status, write_value] = obj.InsertKey(section_name, count_keys+1, key_name, value);
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function [status, write_value] = InsertKey(obj, section_name, key_pos, key_name, value)
            %InsertKey - insert key into the specified position in a given section
            % Using:
            %   [status, write_value] = InsertKey(section_name, key_pos, key_name, value)
            %
            % Output:
            %   status - 1 (true): Success, status - 0 (false): Failed
            %   write_value - 1 (true): Success, status - 0 (false): Failed
            %
            
            if (obj.is_created_configuration)
                if (nargin < 4)
                    error('Not enough input arguments.')
                elseif (nargin == 4)
                    value = '';
                elseif (nargin == 5)
                    if isempty(value)
                        value = '';
                    end
                end
                
                if (~isnumeric(key_pos) && ~isscalar(key_pos))
                    error('Index must be a scalar positive integer')
                else
                    if (key_pos < 1)
                        error('Index must be a positive integer')
                    end
                end
                
                write_value = 0;
                section_name = obj.ValidateSectionName(section_name);
                key_name = obj.ValidateKeyName(key_name);
                sect_num = obj.NameToNumSection(section_name);
                
                if (~isempty(sect_num) && ~isempty(key_name))
                    is_unique_key = obj.IsUniqueKeyName(section_name, key_name);
                    if ~is_unique_key
                        status = false;
                        return;
                    end
                    
                    [keys_indexes, count_keys] = obj.GetKeysIndexes(section_name);
                    if (count_keys > 0)
                        if (key_pos <= count_keys)
                            insert_index = keys_indexes(key_pos);
                        elseif (key_pos > count_keys)
                            insert_index = keys_indexes(end)+1;
                        end
                    else
                        insert_index = obj.indexes_of_sections(sect_num) + 1;
                    end
                    
                    new_data = {key_name, '', ''};
                    
                    obj.config_data_array = insertcell(obj.config_data_array, ...
                        insert_index, new_data);
                    
                    obj.UpdateCountStrings();
                    obj.UpdateSectionsInfo();
                    obj.UpdateEmptyStringsInfo();
                    obj.UpdateCountKeysInfo();
                    
                    if ~isempty(value)
                        write_value = obj.SetValues(section_name, key_name, value);
                    end
                    
                    status = true;
                else
                    status = false;
                end
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function status = RenameKey(obj, section_name, key_name, new_key_name)
            %RenameKey - rename the key in a given section
            % Using:
            %   status = RenameKey(section_name, key_name, new_key_name)       
            %       status - 1 (true): Success
            %       status - 0 (false): Failed
            %
            
            if (obj.is_created_configuration)
                if (nargin < 4)
                    error('Not enough input arguments.')
                end
                
                section_name = obj.ValidateSectionName(section_name);
                key_name = obj.ValidateKeyName(key_name);
                new_key_name = obj.ValidateKeyName(new_key_name);
                sect_num = obj.NameToNumSection(section_name);
                [keys, count_keys] = obj.GetKeys(section_name);
                
                if (~isempty(sect_num) && ~isempty(key_name) && ~isempty(new_key_name) && count_keys > 0)
                    is_unique_key = obj.IsUniqueKeyName(section_name, key_name);
                    if is_unique_key
                        status = false;
                        return;
                    end
                    
                    tf = find(strcmp(key_name, keys), 1, 'last');
                    keys_indexes = obj.GetKeysIndexes(section_name);
                    
                    key_index = keys_indexes(tf);
                    obj.config_data_array{key_index, 1} = new_key_name;
                    
                    status = true;
                else
                    status = false;
                end
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function status = RemoveKey(obj, section_name, key_name)
            %RemoveKey - remove the key from a given section
            % Using:
            %   status = RemoveKey(section_name, key_name)
            %       status - 1 (true): Success
            %       status - 0 (false): Failed
            %
            
            if (obj.is_created_configuration)
                if (nargin < 3)
                    error('Not enough input arguments.')
                end
                
                section_name = obj.ValidateSectionName(section_name);
                key_name = obj.ValidateKeyName(key_name);
                sect_num = obj.NameToNumSection(section_name);
                [keys, count_keys] = obj.GetKeys(section_name);
                
                if (~isempty(sect_num) && ~isempty(key_name) && count_keys > 0)
                    is_unique_key = obj.IsUniqueKeyName(section_name, key_name);
                    if is_unique_key
                        status = false;
                        return;
                    end
                    
                    tf = find(strcmp(key_name, keys), 1, 'last');
                    keys_indexes = obj.GetKeysIndexes(section_name);
                    
                    key_index = keys_indexes(tf);
                    obj.config_data_array(key_index, :) = [];
                    
                    obj.UpdateCountStrings();
                    obj.UpdateSectionsInfo();
                    obj.UpdateEmptyStringsInfo();
                    obj.UpdateCountKeysInfo();
                    
                    status = true;
                else
                    status = false;
                end
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
        %------------------------------------------------------------------
        function status = WriteFile(obj, filename)
            %WriteFile - write to the configuration INI file on disk
            % Using:
            %   status = WriteFile(filename)
            %       status - 1 (true): Success
            %       status - 0 (false): Failed
            %
            
            if (obj.is_created_configuration)
                if (nargin < 1)
                    error('Not enough input arguments.')
                else
                    if ~ischar(filename)
                        error('Requires string input.')
                    else
                        fid = fopen(filename, 'w');
                        
                        if (fid ~= -1)
                            count_str = obj.count_strings;
                            indexes_of_sect = obj.indexes_of_sections;
                            config_data = obj.config_data_array;
                            
                            for k = 1:count_str
                                if isempty(config_data{k,1})
                                    if isempty(config_data{k,3})
                                        fprintf(fid, '\n');
                                    else
                                        comment_str = config_data{k,3};
                                        fprintf(fid, '%s\n', comment_str);
                                    end
                                    
                                elseif ~isempty(indexes_of_sect(indexes_of_sect == k))
                                    if isempty(config_data{k,3})
                                        section_str = config_data{k,1};
                                        
                                        fprintf(fid, '%s\n', section_str);
                                    else
                                        section_str = config_data{k,1};
                                        comment_str = config_data{k,3};
                                        
                                        fprintf(fid, '%s    %s\n', ...
                                            section_str, comment_str);
                                    end
                                    
                                elseif ~isempty(config_data{k,1}) && ...
                                        isempty(indexes_of_sect(indexes_of_sect == k))
                                    if isempty(config_data{k,3})
                                        key_str = config_data{k,1};
                                        val_str = config_data{k,2};
                                        
                                        fprintf(fid, '%s = %s\n', ...
                                            key_str, val_str);
                                    
                                    else
                                        key_str = config_data{k,1};
                                        val_str = config_data{k,2};
                                        comment_str = config_data{k,3};
                                        
                                        fprintf(fid, '%s = %s    %s\n', ...
                                            key_str, val_str, comment_str);
                                    end
                                end
                            end
                        
                            fclose(fid);
                            status = true;
                        else
                            status = false;
                            return;
                        end
                    end
                end
            else
                error('Configuration is not created or not loaded.')
            end
        end
        
    end % public methods
    
    
    methods (Access = 'private')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Private Methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %------------------------------------------------------------------
        function num = NameToNumSection(obj, name)
            %NameToNumSection - get number of section
            
            sections = obj.GetSections();
            indexes = find(strcmp(name, sections));
            
            if ~isempty(indexes)
                % If the section is not unique, then choose the latest
                num = indexes(end);
            else
                num = [];
            end
        end
        
        %------------------------------------------------------------------
        function [keys_indexes, count_keys] = GetKeysIndexes(obj, section_name)
            %GetKeysIndexes - get keys indices from given section
            
            sect_num = obj.NameToNumSection(section_name);
                
            if isempty(sect_num)
                keys_indexes = [];
                count_keys = 0;
                return;
                
            elseif (sect_num == obj.count_sections)
                keys_indexes = ...
                    obj.indexes_of_sections(sect_num)+1:obj.count_strings;
            else
                keys_indexes = ...
                    obj.indexes_of_sections(sect_num)+1:obj.indexes_of_sections(sect_num+1)-1;
            end
            
            indexes_of_empty = obj.indexes_of_empty_str;
            [i_str, empty_indexes] = intersect(keys_indexes, indexes_of_empty);
            
            keys_indexes(empty_indexes) = [];
            keys_indexes = keys_indexes(:);
            count_keys = length(keys_indexes);
        end
        
        %------------------------------------------------------------------
        function UpdateSectionsInfo(obj)
            %UpdateSectionsInfo
            
            keys_data = obj.config_data_array(:,1);
            sect_indexes_cell = regexp(keys_data, '^\[.*\]$');
            obj.indexes_of_sections = find(cellfun(@(x) ~isempty(x), sect_indexes_cell));
            obj.count_sections = length(obj.indexes_of_sections);
        end
        
        %------------------------------------------------------------------
        function UpdateCountKeysInfo(obj)
            %UpdateCountKeys
            
            obj.full_count_keys = obj.count_strings - obj.count_sections - obj.count_empty_strings;
        end
        
        %------------------------------------------------------------------
        function UpdateEmptyStringsInfo(obj)
            %UpdateEmptyStringsInfo
            
            keys_data = obj.config_data_array(:,1);
            indexes_of_empty_cell = strcmp('', keys_data);
            obj.indexes_of_empty_str = find(indexes_of_empty_cell);
            obj.count_empty_strings = length(obj.indexes_of_empty_str);
        end
        
        %------------------------------------------------------------------
        function UpdateCountStrings(obj)
            %UpdateCountStrings
            
            obj.count_strings = size(obj.config_data_array, 1);
        end
        
        %------------------------------------------------------------------
        function status = IsUniqueKeyName(obj, section_name, key_name)
            %IsUniqueKeyName - check whether the name of the key unique
            
            keys = obj.GetKeys(section_name);
            status = ~any(strcmp(key_name, keys));
        end
        
        %------------------------------------------------------------------
        function status = IsUniqueSectionName(obj, section_name)
            %IsUniqueKeyName - check whether the name of the section a unique
            
            sections = obj.GetSections();
            status = ~any(strcmp(section_name, sections));
        end
        
        %------------------------------------------------------------------
        function section_name = ValidateSectionName(obj, section_name)
            %ValidateSectionName - check the name of the section

            if (~ischar(section_name) || size(section_name, 1) > 1)
                error('Requires string input for name of section.')
            else
                section_name = strtrim(section_name);

                if isempty(section_name)
                    section_name = [];
                else
                    sect_indexes_cell = regexp(section_name, '^\[.*\]$', 'once');
                    indexes_cell_comment = regexp(section_name, obj.comment_style, 'once');

                    if ~isempty(indexes_cell_comment)
                        section_name = [];
                        return;
                    end

                    if isempty(sect_indexes_cell)
                        section_name = ['[', section_name, ']'];
                    end
                end
            end
        end

        %------------------------------------------------------------------
        function key_name = ValidateKeyName(obj, key_name)
            %ValidateKeyName - check the name of the key

            if (~ischar(key_name) || size(key_name, 1) > 1)
                error('Requires string input for key name.')
            else
                key_name = strtrim(key_name);
                indexes_cell = regexp(key_name, '^\[.*\]$', 'once');
                indexes_cell_comment = regexp(key_name, obj.comment_style, 'once');

                if isempty(key_name) || ~isempty(indexes_cell) || ~isempty(indexes_cell_comment)
                    key_name = [];
                end
            end
        end
        
    end % private methods
    
end % classdef IniConfig
%--------------------------------------------------------------------------



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tools Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
function config_data = IniFileParser(file_data, comment_style)
    %IniFileParser - parser data from the INI file
    
    % Select the comment in a separate array
    pat = sprintf('^[^%s]+', comment_style);
    comment_data = regexprep(file_data, pat, '');
    
    % Deleting comments
    pat = sprintf('%s.*+$', comment_style);
    file_data = regexprep(file_data, pat, '');
    
    % Select the key value in a separate array
    values_data = regexprep(file_data, '^.[^=]*.', '');
    
    % Select the names of the sections and keys in a separate array
    keys_data = regexprep(file_data, '=.*$', '');
    
    config_data = cell(size(file_data, 1), 3);
    config_data(:,1) = keys_data;
    config_data(:,2) = values_data;
    config_data(:,3) = comment_data;
    
    config_data = strtrim(config_data);
end
%--------------------------------------------------------------------------

%==========================================================================
function values = ParserValues(values)
    %ParserValues - classify the data types and convert them
    
    start_idx = regexp(values, '^[\-\d\s\,\.\:truefalse]+$');
    number_indexes = cellfun(@(x) ~isempty(x), start_idx);
    
    num_values_strs = values(number_indexes);
    num_values = cellfun(@(x) str2num(x), num_values_strs, 'UniformOutput', 0);
    
    empty_num_indexes = cellfun(@(x) isempty(x), num_values);
    num_values(empty_num_indexes) = num_values_strs(empty_num_indexes);
    values(number_indexes) = num_values;
    
    empty_indexes = cellfun(@(x) isempty(x), values);
    values(empty_indexes) = {[]};
end
%--------------------------------------------------------------------------

%==========================================================================
function values = CorrectNumArrayStrings(values)
    %CorrectNumStrings
    
    values = regexprep(values, '^\[', '');
    values = regexprep(values, '\]$', '');
    values = regexprep(values, '\s', ', ');
end
%--------------------------------------------------------------------------

%==========================================================================
function comment_style = ValidateCommentStyle(comment_style)
    %ValidateCommentStyle
    
    if ~ischar(comment_style)
        error('Requires char input for comment style.')
    else
        if (length(comment_style) ~= 1)
            error('The style comments should contain the single character')
        end
    end
end
%--------------------------------------------------------------------------

%==========================================================================
function [keys, values] = ValidateInputsForGetAndSetValues(section, keys, values)
    %ValidateInputsForGetAndSetValues
    
    if ~ischar(section)
        error('First input argument must be a string.')
    end
    
    if ~ischar(keys) && ~iscellstr(keys)
        error('Second input argument must be a string or cell array of strings_indexes.')
    else
        if ischar(keys)
            keys = {keys};
        end
    end
    
    if ~iscell(values)
        values = {values};
    end
    
    if (numel(keys) ~= numel(values))
        error('Number of elements in the second and third input argument must be equal')
    end
    
    keys = keys(:);
    values = values(:);
end
%--------------------------------------------------------------------------

%==========================================================================
function B = insertcell(C, index, data)
    %insertcell - insert a new cell or several cells in a two-dimensional
    % array of cells on the index
    
    M = index;
    N = 1;
    
    if ~iscell(data)
        data = {data};
    end
	
    [m, n] = size(C);
    [md, nd] = size(data);
    
    mm = m + md;
    nt = nd + N - n - 1;
    
    if (nt > 0)
        nn = n + nt;
    else
        nn = n;
    end
    
    B = cell(mm, nn);
    B(M:M+md-1, N:N+nd-1) = data;
    
    ia = 1:M-1;
    ib = M+md:mm;
    j = 1:n;
    
    ibb = [ia, ib];
    if (m > 2)
        t = 0;
    else
        t = 1;
    end
    icc = [ia, M:length(C)-t];
    B(ibb, j) = C(icc, j);
end
%--------------------------------------------------------------------------

