function save_adress_data(dataStr, SaveFullName, isoslice_args)
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
if nargin < 3; isoslice_args = []; end
if isempty(isoslice_args); isoslice_args = []; end
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

%% 3 - Saving a text file with all the info
if isfield(dataStr, 'meta')
    fileID = fopen(char(save_fullfile(1:end-4)+"_info.txt"), 'w');
    fprintf(fileID, dataStr.meta.info);
    fclose(fileID);
end

%% 4 - Saving a figure snapshot with the data
% (A) -- If the data to be saved is an IsoSlice / IsoCut /IsoScan form
if isfield(dataStr, 'IsoType')
    pp = plot_props();
% For a IsoSlice data structure
    if dataStr.IsoType == "IsoE" || dataStr.IsoType == "IsoK"
        fig = figure(); fig.Position(3) = pp.fig5x4(1); fig.Position(4) = pp.fig5x4(2);
        ImData(dataStr.XSlice, dataStr.YSlice, dataStr.DSlice);
        % - Formatting the figure
        axis([min(dataStr.XSlice(:)), max(dataStr.XSlice(:)), min(dataStr.YSlice(:)), max(dataStr.YSlice(:))]);
        minC = min(dataStr.DSlice(:)); maxC = max(dataStr.DSlice(:));
        caxis([minC, maxC]);
        img_props(); colbar_props();
         % - Re-labelling the axes depending on what slice is taken
        if dataStr.IsoType== "IsoE"
            if dataStr.Type == "Eb(kx,ky)"
                if isfield(dataStr, 'kx'); xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  \tau (^{\circ}) $$', 'Interpreter', 'latex');
                end
            elseif dataStr.Type == "Eb(kx,kz)"
                if isfield(dataStr, 'kx'); xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
                end
            elseif dataStr.Type == "Eb(k,i)"
                if isfield(dataStr, 'kx'); xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); ylabel('$$ \bf  Index $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  Index $$', 'Interpreter', 'latex'); ylabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
                end
            end
        elseif dataStr.IsoType== "IsoK"
            ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex'); 
            if dataStr.Type == "Eb(kx,ky)"
                if isfield(dataStr, 'kx'); xlabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex');
                else; xlabel('$$ \bf  \tau (^{\circ}) $$', 'Interpreter', 'latex');
                end
            elseif dataStr.Type == "Eb(kx,kz)"
                if isfield(dataStr, 'kx'); xlabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  hv (eV) $$', 'Interpreter', 'latex');
                end
            elseif dataStr.Type == "Eb(k,i)"
                if isfield(dataStr, 'kx'); xlabel('$$ \bf  Index $$', 'Interpreter', 'latex'); 
                else; xlabel('$$ \bf  Index $$', 'Interpreter', 'latex');
                end
            end
        end
        title(sprintf(string(dataStr.FileName) + "; %s; [%.2f,%.2f]", dataStr.IsoType, dataStr.IsoWin(1),dataStr.IsoWin(2)), 'interpreter', 'none', 'fontsize', 8);
        print(fig, char(save_fullfile(1:end-4)+"_snap.png"), '-dpng');

% For a IsoScan data structure
    elseif dataStr.IsoType == "Scan"
        fig = figure(); fig.Position(3) = pp.fig5x4(1); fig.Position(4) = pp.fig5x4(2);
        ImData(dataStr.XScan, dataStr.YScan, dataStr.DScan);
        img_props([], string(dataStr.xField));
        colbar_props();
        title(sprintf(string(dataStr.FileName) + "; scan %i", dataStr.IsoWin), 'interpreter', 'none', 'fontsize', 8);
        print(fig, char(save_fullfile(1:end-4)+"_snap.png"), '-dpng');

% For a IsoCut data structure
    elseif dataStr.IsoType == "EDC" || dataStr.IsoType == "MDC"
        fig = figure(); hold on;
        if dataStr.IsoType == "MDC"
            fig.Position(3) = pp.fig16x9(1); 
            fig.Position(4) = pp.fig16x9(2);
            plot(dataStr.XCut, dataStr.DCut, 'k.-', 'color', [0 1 0], 'LineWidth', 1.5);
            gca_props(); box off;
            % Axis labels and limits
            if dataStr.xField == "raw_tht" || dataStr.xField == "tht"
                xlabel('$$ \bf  \theta (^{\circ}) $$', 'Interpreter', 'latex');
            elseif dataStr.xField == "kx"
                xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
            end
            ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
            % - Plotting the x- and y-axes
            xl = [min(dataStr.XCut(:)), max(dataStr.XCut(:))]; yl = ylim;
            axis([xl(1), xl(2), yl(1), 1.05*yl(2)]);
        elseif dataStr.IsoType== "EDC"
            fig.Position(3) = pp.fig16x9(2); 
            fig.Position(4) = pp.fig16x9(1);
            plot(dataStr.DCut, dataStr.XCut, 'k.-', 'color', [0 0 1], 'LineWidth', 1.5);
            gca_props(); box off;
            % - Axis labels and limits
            ylabel('$$ \bf  E_{B} (eV) $$', 'Interpreter', 'latex');
            xlabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
            % - Plotting the x- and y-axes
            xl = xlim; yl = [min(dataStr.XCut(:)), max(dataStr.XCut(:))];
            axis([xl(1), 1.05*xl(2), yl(1), yl(2)]);
        end
        title(sprintf(string(save_filename) + "; %s; [%.2f,%.2f]", dataStr.IsoType, dataStr.IsoWin(1),dataStr.IsoWin(2)), 'interpreter', 'none', 'fontsize', 8);
        print(fig, char(save_fullfile(1:end-4)+"_snap.png"), '-dpng');
    % For a IsoCut data structure
    elseif dataStr.IsoType == "XPS"
        fig = view_xps(dataStr);
        print(fig, char(save_fullfile(1:end-4)+"_snap.png"), '-dpng');
        % -- If fit is available
        if isfield(dataStr, 'fit')
            fig = view_xps_fit(dataStr.fit); dataStr.fit
            print(fig, char(save_fullfile(1:end-4)+"_fit.png"), '-dpng');
        end
    end
    
% (B) -- Else, for the usual ARPES data structure
elseif isfield(dataStr, 'Type')
    if dataStr.Type == "Eb(k)" || dataStr.Type == "Eb(kx,ky)" || dataStr.Type == "Eb(kx,kz)" || dataStr.Type == "Eb(k,i)"
        fig = view_data(dataStr, isoslice_args);
        print(fig, char(save_fullfile(1:end-4)+"_snap.png"), '-dpng');
    elseif dataStr.Type == "arpes2model fit"
        % -- If fit is available
        if isfield(dataStr, 'fit')
            fig = view_arpes2model_fit(dataStr.fit); dataStr.fit
            print(fig, char(save_fullfile(1:end-4)+"_fit.png"), '-dpng');
        end
    end
end

%% Close wait-bar
waitbar(1,'Save complete!'); 
end