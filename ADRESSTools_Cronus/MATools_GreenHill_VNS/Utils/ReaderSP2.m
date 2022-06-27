function [data,comments,transform] = ReaderSP2(filename)
    data = struct('image',[], ...
        'raw', [], ...   % raw image without correction
        'erange',[], ... % in the format of [EkinStart, EkinEnd]
        'arange',[], ... % in the format of [aStart, aEnd]
        'edelta',[], ... % in the format of [EkinStart, EkinDelta]
        'adelta',[]  ... % in the format of [aStart, aDelta]
        );
    d = readsp2(filename);
    if (isempty(d.image))
        return
    end
    data.image = d.image;
    data.raw = d.raw;
    comments = d.comment;
    
    % Parse ERange and aRange
    ini = IniConfig();
    ini.ReadStream(comments,'#');
    
    value = ini.GetValues('Transform','ERange','');
    erange = value{1};
    delta = (erange(2)-erange(1))/(size(data.image,2)-1);
    data.erange = erange;
    data.edelta = [erange(1) delta];

    value = ini.GetValues('Transform','aRange','');
    arange = value{1};
    delta = (arange(2)-arange(1))/(size(data.image,1)-1);
    data.arange = arange;
    data.adelta = [arange(1) delta];
    
    % Parse Transform section
    transform = [];
    keys = {'lensmode', 'Ek', 'Ep', 'WF', 'ERange', 'aRange', 'De1', 'aInner', ...
            'eShift', 'Da1', 'Da3', 'Da5', 'Da7', 'pixel_size', 'magnification', ...
            'offset', 'interpolation_mode', 'intensity_scaling', 'aUnit', ...
            'eGrid', 'aGrid', 'source'};
    for key=keys
        keyname = key{1};
        value = ini.GetValues('Transform', keyname, '');
        value = value{1};
        transform.(keyname) = value;
    end
end
