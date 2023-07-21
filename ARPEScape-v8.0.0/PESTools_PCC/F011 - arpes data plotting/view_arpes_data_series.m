function view_arpes_data_series(dataStr, caxLims)
% view_arpes_data_series(dataStr, caxLims)
%	General function to view a series of all scans from 3D ARPES data. Make
%   sure ONLY Eb(kx,ky) or Eb(kx,kz) data is used here. This function plots 
%   the ARPES data in the form of D(X,Y) for all the scan parameter values 
%   (tltM or hv) in the form of a  downsampled image series. This function 
%   ensures that the most recent processing is shown in the figure.
%
%   REQ. FUNCTIONS:
%   -   [xField, yField, zField, dField] = find_data_fields(dataStr);
%   -   [h=] ImData(X,Y,Z[,style])
%
%   IN:
%   -   dataStr:          data structure of the ARPES data.
%   -   caxLims:          1x2 vector in normalised units of the colorscale of each image. 
%
%   OUT: (none)

%% Default parameters
if nargin < 2; caxLims = [0, 1]; end
if isempty(caxLims); caxLims = [0, 1]; end

disp('-> view series of all scans...')
wbar = waitbar(0., 'Plotting a view series of all scans...', 'Name', 'view_data_series');

%% - 1 - Initialising plot parameters
pp = plot_props();
% - Extracting the fields to be used with most recent processing
[xField, yField, zField, dField] = find_data_fields(dataStr);
% - Each series is an nRows x nCols sub-plot
nRows = 3; nCols = 3;

%% - 2 - Down-sampling the form of the data
step_size = 2;
zMidIndx = ceil(size(dataStr.(dField), 3)/2);

%% - 3 - Plotting a series of Eb(k) Images
% Filing through all the pages of the figures
for ipage=1:ceil(size(dataStr.(dField), 3)/nRows/nCols)
    figure_cell{ipage} = figure(); orient tall; 
    set(figure_cell{ipage}, 'Position', [80, 80, pp.fig8x8(1), pp.fig8x8(2)], 'Name', dataStr.H5file);
    % Filing through all the frames in a single figure
    for iframe=1:nRows*nCols
        % - Finding the n'th data object to be plotted
        n = (ipage-1)*nRows*nCols+iframe;
        waitbar(n/size(dataStr.(dField), 3), wbar, 'Plotting a view series of all scans...', 'Name', 'view_data_series');
        % - If n > number of plotting elements, break the statement
        if n > size(dataStr.(dField), 3)
            break; 
        else
            % -- Plotting the n'th data object
            subplot(nRows,nCols,iframe); 
            img_props([], string(xField));
            xlabel(''); ylabel('');
            % -- The plot depends on how far in the analysis 
            % - 2.1 - EbAlign->Normalise->kConvert fields
            if isfield(dataStr, 'kx')
                ImData(dataStr.(xField)(1:step_size:end,1:step_size:end,n), dataStr.(yField)(1:step_size:end,1:step_size:end,n), dataStr.(dField)(1:step_size:end,1:step_size:end,n));
                title(sprintf('scan = %.2f', dataStr.(zField)(zMidIndx,zMidIndx,n)), 'fontsize', 10);
            % - 2.2 - EbAlign->Normalise fields
            elseif isfield(dataStr, 'data')
                ImData(dataStr.(xField)(1:step_size:end,1:step_size:end,n), dataStr.(yField)(1:step_size:end,1:step_size:end,n), dataStr.(dField)(1:step_size:end,1:step_size:end,n));
                title(sprintf('scan = %.2f', dataStr.(zField)(n)), 'fontsize', 10);
            % - 2.3 - EbAlign fields
            elseif isfield(dataStr, 'eb')
                ImData(dataStr.(xField)(1:step_size:end,1:step_size:end,n), dataStr.(yField)(1:step_size:end,1:step_size:end,n), dataStr.(dField)(1:step_size:end,1:step_size:end,n));
                title(sprintf('scan = %.2f', dataStr.(zField)(n)), 'fontsize', 10);
            % - 2.4 - Raw, unprocessed data fields
            else
                ImData(dataStr.(xField)(1:step_size:end), dataStr.(yField)(1:step_size:end), dataStr.(dField)(1:step_size:end,1:step_size:end,n));
                title(sprintf('scan = %.2f', dataStr.(zField)(n)), 'fontsize', 10);
            end
            % - 2.5 - Formatting the figure
            minC = min(min(min(dataStr.(dField)(1:step_size:end,1:step_size:end,iframe))));
            maxC = max(max(max(dataStr.(dField)(1:step_size:end,1:step_size:end,iframe))));
            caxis([caxLims(1), caxLims(2)].*maxC);
            
            % -- Remove the colorbar for all plots, except the first one
            if iframe == 1; colbar_props([0.925 0.80 0.020 0.125]); end
        end
    end
end
close(wbar);
end