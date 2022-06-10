function [cubeStr, fig] = extract_isoCube(dataStr, isocube_args, plot_results)
% [cubeStr, fig] = extract_isoCube(dataStr, isocube_args, plot_results)
%   This function plots a 3D ARPES data data cube in the form
%   of D(X, Y, Z). General function to view 3D ARPES data via a data cube 
%   approach
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   [h=] ImData(X,Y,Z[,style])
%
%   IN:
%   -   dataStr:       	data structure of the ARPES data.
%   -   isocube_args:  	1x2 cell of {cLimMin, smoothVal}.
%   -   plot_results: 	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -   cubeStr:    	new MATLAB data-structure.
%   -   fig:            figure output.


%% Default parameters
if nargin < 3; plot_results = 1; end
if nargin < 2; isocube_args = cell(1,2); plot_results = 1; end
if isempty(isocube_args); isocube_args = cell(1,2); end
if isempty(plot_results); plot_results = 1; end
pp = plot_props();
fig = [];
% - Extracting the fields to be used with most recent processing
[xField, yField, zField, dField] = find_data_fields(dataStr);

%% - 1 - Initialising input parameters
cLimMin         = isocube_args{1}; if isempty(cLimMin);     cLimMin = 0.25; end     % minimum threshold, below which the data points are transparent
smoothVal       = isocube_args{2}; if isempty(smoothVal);   smoothVal = 5;  end     % FWHM of a Gaussian smoothing operation on the 3D cube

%% 1 - Initialising and permuting the variable / data arrays to be consistent in 3D
% - 1.1 - EbAlign->Normalise->kConvert fields
if isfield(dataStr, 'kx')
    v = permute(dataStr.(dField), [3 2 1]);
    x = permute(dataStr.(xField), [3 2 1]);
    y = permute(dataStr.(yField), [3 2 1]);
    z = permute(dataStr.(zField), [3 2 1]);
% - 1.2 - EbAlign->Normalise fields
elseif isfield(dataStr, 'data')
    v = permute(dataStr.(dField), [3 2 1]);
    x = permute(dataStr.(xField), [3 2 1]);
    y = permute(dataStr.(yField), [3 2 1]);
    % - Converting the scan variable into 3D
    z = repmat(dataStr.(zField), [size(dataStr.(dField),1),1,size(dataStr.(dField),2)]);
    z = permute(z, [2 3 1]);
% - 1.3 - EbAlign fields
elseif isfield(dataStr, 'eb')
    v = permute(dataStr.(dField), [3 2 1]);
    x = permute(dataStr.(xField), [3 2 1]);
    y = permute(dataStr.(yField), [3 2 1]);
    % - Converting the scan variable into 3D
    z = repmat(dataStr.(zField), [size(dataStr.(dField),1),1,size(dataStr.(dField),2)]);
    z = permute(z, [2 3 1]);
% - 1.4 - Raw, unprocessed data fields
else
    v = permute(dataStr.(dField), [3 2 1]);
    % - Converting the theta variable into 3D
    x = repmat(dataStr.(xField), [size(dataStr.(dField),1),1,size(dataStr.(dField),3)]);
    x = permute(x, [3 2 1]);
    % - Converting the eb variable into 3D
    y = repmat(dataStr.(yField), [1,size(dataStr.(dField),2),size(dataStr.(dField),3)]);
    y = permute(y, [3 2 1]);
    % - Converting the scan variable into 3D
    z = repmat(dataStr.(zField), [size(dataStr.(dField),1),1,size(dataStr.(dField),2)]);
    z = permute(z, [2 3 1]);
end

%% 2 - Interpolating the data with a Gaussian if necessary
% - If necessary, smooth the data beforehand
v = smooth3(v, 'gaussian', smoothVal);

%% 2.0 - Initialising the figure
if plot_results == 1
    fig = figure(); hold on;
    fig.Position(1) = pp.figpos(1);
    fig.Position(2) = pp.figpos(2);
    fig.Position(3) = pp.fig5x4(1); 
    fig.Position(4) = pp.fig5x4(2);
    % - For a 3D Eb(kx,ky) scan type
    if dataStr.Type == "Eb(kx,ky)"
        % -- Plotting the Eb(kx,ky) figure object
        if isfield(dataStr, 'kx'); set(fig, 'Name', 'Iso-Scan - Eb(kx,ky): A, N, K');
        elseif isfield(dataStr, 'data'); set(fig, 'Name', 'Iso-Scan - Eb(kx,ky): A, N');
        elseif isfield(dataStr, 'eb'); set(fig, 'Name', 'Iso-Scan - Eb(kx,ky): A');
        else; set(fig, 'Name', 'Iso-Scan - Eb(kx,ky): Raw');
        end
    % - For a 3D Eb(kx,kz) scan type
    elseif dataStr.Type == "Eb(kx,kz)"
        % -- Plotting the Eb(kx,kz) figure object
        if isfield(dataStr, 'kx'); set(fig, 'Name', 'Iso-Scan - Eb(kx,kz): A, N, K');
        elseif isfield(dataStr, 'data'); set(fig, 'Name', 'Iso-Scan - Eb(k): A, N');
        elseif isfield(dataStr, 'eb'); set(fig, 'Name', 'Iso-Scan - Eb(kx,kz): A');
        else; set(fig, 'Name', 'Iso-Scan - Eb(kx,kz): Raw');
        end
    end
    % - Plotting the iso surfaces of the cube
    p1L = patch(isocaps(z, x, y, v, cLimMin*max(v(:))), 'FaceColor', 'interp', 'EdgeColor', 'none');
    p2L = patch(isosurface(z, x, y, v, cLimMin*max(v(:))),'FaceColor', 'black', 'EdgeColor', 'none');
    isonormals(v, p2L);
    % - Formatting the figure
    img_props([], string(xField));
    colbar_props([0.925 0.80 0.020 0.125]);
    % - View orientation
    view([-65, 12]);
    axis([min(z(:)), max(z(:)), min(x(:)), max(x(:)), min(y(:)), max(y(:))]);
    % - Defining the axes properties
    ax = gca;
    ax.Color        = [1, 1, 1];
    ax.XColor       = [0 0 0]; 
    ax.YColor       = [0 0 0]; 
    ax.ZColor       = [0 0 0];
    ax.Box          = 'on'; 
    ax.BoxStyle     = 'full'; 
    ax.Layer        = 'Top';
    % - Axis labels and limits
    if dataStr.Type == "Eb(kx,ky)" 
        if isfield(dataStr, 'kx'); ylabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); xlabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex');
        else; ylabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex'); xlabel('$$ \bf  \tau (^{\circ}) $$', 'Interpreter', 'latex');
        end
    elseif dataStr.Type == "Eb(kx,kz)" 
        if isfield(dataStr, 'kx'); ylabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); xlabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex');
        else; ylabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex'); xlabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
        end
    end
    zlabel('$$ \bf  E_B (eV) $$', 'Interpreter', 'latex');
end

%% 6 - Appending data to MATLAB data structure
cubeStr                 = struct();
cubeStr.isocube_args   	= isocube_args;
cubeStr.xField          = string(xField);
cubeStr.FileName        = dataStr.FileName;
cubeStr.Type            = dataStr.Type;
cubeStr.IsoType         = "Cube";
cubeStr.XCube         	= x;
cubeStr.YCube         	= y;
cubeStr.ZCube         	= z;
cubeStr.DCube          	= v;
end