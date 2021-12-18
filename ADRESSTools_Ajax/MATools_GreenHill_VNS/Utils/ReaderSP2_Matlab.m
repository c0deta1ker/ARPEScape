function [data,comments] = ReaderSP2(filename)
    % Read SP2 files.
    %   It read the image as well as the scale information.
    
    data = struct('image',[], ... 
        'erange',[], ... % in the format of [EkinStart, EkinEnd]
        'arange',[], ... % in the format of [aStart, aEnd]
        'edelta',[], ... % in the format of [EkinStart, EkinDelta]
        'adelta',[]  ... % in the format of [aStart, aDelta]
        );


    % Read image data
    % This is much slower than read with textscan
%     try
%         data.image = imread(filename);
%         lastwarn('');
%     catch ME
%         disp(ME.message)
%         return
%     end
    
    % Valid SP2 file?
    file = fopen(filename);
    line = fgetl(file);
    if ~strcmp(line, 'P2')
        fclose(file);
        return
    end

    % Read comments
    comments = {};
    line = fgetl(file);
    while ischar(line)
        if line(1)~='#', break; end
        comments = vertcat(comments,regexprep(line, '^#(.*)', '$1'));
        line = fgetl(file);
    end

    % Read dimension information
    dims = sscanf(line, '%d %d %d');
    % Read data
    raw = textscan(file, '%d');
    data.image = double(reshape(raw{1}, dims(1), dims(2))');
    % Close file
    fclose(file);
    
    % Parse ERange and aRange
    ini = IniConfig();
    ini.ReadStream(comments,'#');

    value = ini.GetValues('Transform','ERange','');
    range = value{1};
    delta = (range(2)-range(1))/(size(data.image,2)-1);
    data.erange = range;
    data.edelta = [range(1) delta];

    value = ini.GetValues('Transform','aRange','');
    range = value{1};
    delta = (range(2)-range(1))/(size(data.image,1)-1);
    data.arange = range;
    data.adelta = [range(1) delta];
    
end

