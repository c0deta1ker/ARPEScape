function [data, axes, note] = ReaderHDF5(fname)
    % Read HDF5 file saved by SmartGUI
    %   data is the possible multidimensional matrix.
    %   axes is a cell object. Each row holds the scale array for
    %   corresponding dimension.
    axes = {};
    % data dimensions are reversed
    data  = h5read(fname, '/Matrix');
    data  = permute(data, ndims(data):-1:1);
    % notes
    note  = h5readatt(fname, '/Matrix', 'IGORWaveNote');
    % scale information for each dimension
    scale = h5readatt(fname, '/Matrix', 'IGORWaveScaling');
    [~, ncols] = size(scale);
    for n=2:ncols
        axis = scale(2,n)+[0:size(data,n-1)-1]*scale(1, n);
        axes = vertcat(axes, axis);
    end
end