%% *Initialising the script*
close all; clear all;
% Sample & Dataset properties
filename            = "id20230805_ExXPS";
beamtime_id         = "20230805";
dataset_id          = "Example ARPES data";
sample_id           = "('London')";
sample_stack_info   = "";
measurement_notes   = "Sample measured @ADRESS, 77 K.";
% Defining the paths of the model, data and save directories
path_data   = 'D:\OneDrive\PCC_MS_Project\00_SLS_Software\ARPEScape\ARPEScape-v8.1.0\PESTools_PCC\Examples\Example_Data\ADRESS_data\';
path_save   = 'D:\OneDrive\PCC_MS_Project\00_SLS_Software\ARPEScape\ARPEScape-v8.1.0\PESTools_PCC\Examples\Example_Data_Saved\';

%% *(1) XPS - Data Processing*
%% (1.1) - Load in data
close all; 
% Load in XPS data
file_names      = {...
    '011_XPS_sAl_350.h5',...
    '012_XPS_sAl_450.h5',...
    '013_XPS_sAl_550.h5',...
    '014_XPS_sAl_650.h5',...
    '015_XPS_sAl_750.h5',...
    '016_XPS_sAl_850.h5',...
    '017_XPS_sAl_950.h5',...
    '018_XPS_sAl_1050.h5',...
    '019_XPS_sAl_1150.h5',...
    };
xps_dat = cell(1,length(file_names));
for i = 1:length(file_names); xps_dat{i} = load_adress_data(file_names{i}, path_data); end
view_arpes_data(xps_dat{1});
% Load in Reference data
file_names_ef      = {...
    '011_XPS_350_EF.h5',...
    '012_XPS_450_EF.h5',...
    '013_XPS_550_EF.h5',...
    '014_XPS_650_EF.h5',...
    '015_XPS_750_EF.h5',...
    '016_XPS_850_EF.h5',...
    '017_XPS_950_EF.h5',...
    '018_XPS_1050_EF.h5',...
    '019_XPS_1150_EF.h5',...
    };
xps_ref = cell(1,length(file_names_ef));
for i = 1:length(file_names_ef); xps_ref{i} = load_adress_data(file_names_ef{i}, path_data); end
view_arpes_data(xps_ref{1});
%% (1.2) - Binding energy alignment to the reference data
close all; 
% - AUTOMATIC COARSE ALIGNMENT TO THE EDGE OF THE REFERENCE
% -- Coarse alignment of the energy
alignType   = "align2ref";
scan_indxs	= [];
eWin        = 1.00;
dEWin       = 4.00;
dESmooth    = [];
feat        = 'edge';
align_args  = {alignType, scan_indxs, eWin, dEWin, dESmooth, feat};
for i = 1:length(xps_dat)
    [xps_dat{i}, xps_ref{i}] = align_energy(xps_dat{i}, align_args, xps_ref{i}, 1);
end
% - FITTING TO THE FDD WITH FINE ALIGNMENT
alignType   = "fit2ef"; 
scan_indxs	= [];
eWin        = 0.00;
dEWin       = 0.50;
dESmooth    = [];
feat        = 'edge';
align_args  = {alignType, scan_indxs, eWin, dEWin, dESmooth, feat};
for i = 1:length(xps_dat)
    xps_dat{i}  = align_energy(xps_dat{i}, align_args, xps_ref{i}, 1); 
end
view_arpes_data(xps_dat{1});
%% (1.3) - Angle-integration of data to extract XPS curve
close all; 
for i = 1:length(xps_dat)
    cutWin      = 0.95*[min(xps_dat{i}.raw_tht(:)), max(xps_dat{i}.raw_tht(:))];      % Integration window of the cut to be made.
    xps_args    = {cutWin};
    [xps_dat{i}, ~] = pes_extract_data(xps_dat{i}, xps_args, 0); 
end
view_pes_data(xps_dat);
%% (1.4) - Extracting peak positions vs measurement number
close all; 
% -- Extracting the energy and intensity of the peak positions
xWin       = -17.5 + 0.5*[-1,1];
for i = 1:length(xps_dat)
    [xps_dat{i}.peak_loc, xps_dat{i}.peak_int] = find_peak_loc(xps_dat{i}.xdat, xps_dat{i}.ydat, xWin, 'spline', 1); 
end
% -- Plotting the comparison
fig1 = figure(); fig1.Position(3) = 2.0*500; fig1.Position(4) = 1.0*300;
% -- Plotting the Peak Eb position vs Photon Energy
subplot(121); hold on;
for i = 1:length(xps_dat); plot(i, xps_dat{i}.peak_loc, 'kx-', 'markersize', 10); end
gca_props(); xlim([0.5, length(xps_dat)+0.5]);
title('Peak EB vs index');
ylabel('Peak $$E_{B}$$ [eV]', 'interpreter', 'latex');
% -- Plotting the Peak Intensity vs Photon Energy
subplot(122); hold on;
for i = 1:length(xps_dat); plot(i, xps_dat{i}.peak_int, 'rx-', 'markersize', 10); end
gca_props(); xlim([0.5, length(xps_dat)+0.5]);
title('Peak Intensity vs index');
ylabel('Peak Intensity [counts]', 'interpreter', 'latex');
%% (1.5) - Normalising the XPS data
close all; 
for i = 1:length(xps_dat)
    eWin        = xps_dat{i}.peak_loc;
    dEWin       = 1.0;
    norm_args   = {eWin, dEWin};
    xps_dat{i}.ydat = pes_norm2peak(xps_dat{i}.xdat, xps_dat{i}.ydat, norm_args, 1);
end
view_pes_data(xps_dat);
%% (1.6) - Cropping the XPS data around ROI
close all; 
ebLims = xps_dat{1}.peak_loc + [-3.5, 3.0];
for i = 1:length(xps_dat)
    [~, xps_dat{i}.xdat, xps_dat{i}.ydat] = data_crop2D([], xps_dat{i}.xdat, xps_dat{i}.ydat, [], ebLims);
end
view_pes_data(xps_dat);
%% (1.7) - Final plot of XPS data
close all; 
view_pes_data(xps_dat); xps_dat{1}

%% *(2) FitCL2PEAK: Extracting the Peak Position of the XPS data*
close all; 
%% (2.1) - Initialising XPS data
xps_fit               = struct();
xps_fit.FileName      = file_names;
xps_fit.repeats       = length(xps_dat);
xps_fit.data          = xps_dat;
%% (2.2) - Peak position determination
xWin       = xps_fit.data{1}.peak_loc + 0.10*[-1,1];
xps_fit.ECL = []; xps_fit.ECL_int = [];
for i = 1:length(xps_dat)
    [xps_fit.ECL(i), xps_fit.ECL_int(i)] = find_peak_loc(xps_dat{i}.xdat, xps_dat{i}.ydat, xWin, 'G', 1); 
end
%% (2.3) - Finding the best estimate of ECL across the dataset
close all; 
ECL     = zeros(1,length(xps_fit.ECL)); 
for i = 1:length(xps_dat)
    ECL(i)      = xps_fit.ECL(i);
end
% -- Isolating any potential outliers
XX              = 1:length(ECL);
TF              = isoutlier(ECL); % TF(1:end) = 0;
XX_adj{1}       = XX(~TF);
XX_adj{2}       = XX(TF);
ECL_adj{1}      = ECL(~TF);
ECL_adj{2}      = ECL(TF);
% -- Plotting the result of outlier analysis
figure(); hold on;
plot(XX_adj{1}, ECL_adj{1}, 'kx');
plot(XX_adj{2}, ECL_adj{2}, 'rx');
% -- Final variable assignments
ECL             = ECL(~TF);
ECL_mu          = mean(ECL)
ECL_3sig        = 3*std(ECL)
%% (2.4) - Saving corresponding .mat file
field_name = (filename+"_FitCL2PEAK");
% Defining the MATLAB data structure
dataStr.(field_name)                        = struct();
% - Measurement / Sample Information
dataStr.(field_name).filename               = field_name;
dataStr.(field_name).fit_method             = "FitCL2PEAK";
dataStr.(field_name).beamtime_id            = beamtime_id;
dataStr.(field_name).dataset_id             = dataset_id;
dataStr.(field_name).sample_id              = sample_id;
dataStr.(field_name).sample_stack_info      = sample_stack_info;
% - Best Fit XPS Results
% -- XPS data and corresponding best fits
dataStr.(field_name).num_of_repeats         = length(xps_fit.data(~TF));
dataStr.(field_name).num_of_outliers        = length(xps_fit.data(TF));
dataStr.(field_name).xps_data               = xps_fit.data(~TF);
% -- ECL estimates & measurement uncertainty
dataStr.(field_name).ECL                    = ECL;
dataStr.(field_name).ECL_mu                 = ECL_mu;
dataStr.(field_name).ECL_3sig               = ECL_3sig;
% - Saving the data
save_adress_data(dataStr.(field_name), path_save + field_name);
%% (2.5) - Saving corresponding .txt file
% Creating an information string
info = "";
info = info + sprintf("Filename: \t\t %s \n",   dataStr.(field_name).filename);
info = info + sprintf("Fit Method: \t %s \n",   dataStr.(field_name).fit_method);
info = info + sprintf("Beamtime ID: \t %s \n",  dataStr.(field_name).beamtime_id);
info = info + sprintf("Dataset ID: \t\t %s \n",  dataStr.(field_name).dataset_id);
info = info + sprintf("Sample ID: \t\t %s \n",  dataStr.(field_name).sample_id);
info = info + sprintf("Sample stack: \t %s \n", dataStr.(field_name).sample_stack_info);
info = info + sprintf("Repeats: \t\t %i \n",    dataStr.(field_name).num_of_repeats);
info = info + sprintf("Outliers: \t\t %i \n",    dataStr.(field_name).num_of_outliers);
info = info + sprintf("ECL: \t\t\t %.3f (%.3f) eV \n", dataStr.(field_name).ECL_mu, dataStr.(field_name).ECL_3sig);
info = info + sprintf("Notes: \t\t %s", measurement_notes);
fprintf(info);
% Writing the string to a text file
fileID = fopen(char(path_save + dataStr.(field_name).filename + ".txt"),'w');
fprintf(fileID,info);
fclose(fileID);
%% (2.6) - Saving corresponding figure summary
pp  = plot_props();
% Initialising the figure
fig = figure();
fig.Position(3) = pp.fig6x6(1);
fig.Position(4) = pp.fig6x6(2);
% Plotting the XPS data
subplot(3,1,[1,2]); hold on;
% - Plotting all of the curve components
if length(dataStr.(field_name).xps_data) == 1; cols = [0 0 0];
else; cols = jet(length(dataStr.(field_name).xps_data));
end
x_vals = {}; y_vals = {};
for i = 1:length(dataStr.(field_name).xps_data)
    plot(dataStr.(field_name).xps_data{i}.xdat, dataStr.(field_name).xps_data{i}.ydat,...
        'k-', 'color', cols(i,:), 'linewidth', 1.5);
    x_vals{i} = dataStr.(field_name).xps_data{i}.xdat;
    y_vals{i} = dataStr.(field_name).xps_data{i}.ydat;
end
% - Plotting best fit and 3sig error
patch(dataStr.(field_name).ECL_mu + [-1, -1, 1, 1, -1]*dataStr.(field_name).ECL_3sig,...
    [0, 1, 1, 0, 0],...
    [0 0 0], 'FaceColor', 'k', 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
line([1,1]*dataStr.(field_name).ECL_mu, [0,1], 'Color', 'k', 'LineWidth', 1, 'Linestyle', '-');
% - Formatting the figure
gca_props(); 
xlabel('$$ \bf E_B\ -\ E_F\ (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf Intensity\ (arb.) $$', 'Interpreter', 'latex');
axis([min(cell2mat(x_vals(:))), max(cell2mat(x_vals(:))), min(cell2mat(y_vals(:))), 1.1*max(cell2mat(y_vals(:)))]);
title_text = strrep(dataStr.(field_name).filename, "_", "-");
title(title_text, 'interpreter', 'none');
lgnd = {}; for i = 1:length(dataStr.(field_name).xps_data); lgnd{end+1} = string(dataStr.(field_name).xps_data{i}.hv) + " eV"; end
legend(lgnd,'Location','best');
% Plotting the ECL versus measurement number
subplot(3,1,3); hold on;
XX = 1:length(dataStr.(field_name).ECL);
for i = 1:length(dataStr.(field_name).xps_data)
    plot(XX(i), dataStr.(field_name).ECL(i), 'ko',...
    'color', 'k', 'markerfacecolor', cols(i,:), 'linewidth', 0.5,'markersize', 8);
end
errorbar(mean(XX(:)), dataStr.(field_name).ECL_mu, dataStr.(field_name).ECL_3sig, dataStr.(field_name).ECL_3sig, 0, 0, 'ks-',...
    'color', [0 0 0], 'markerfacecolor', pp.col.fit{2}, 'markersize', 10, 'linewidth', 1.5);
% - Formatting the figure
gca_props(); % grid on; 
axis([...
    0.75, max(XX)+0.25,...
    dataStr.(field_name).ECL_mu - 1.25*dataStr.(field_name).ECL_3sig-0.01,...
    dataStr.(field_name).ECL_mu + 1.25*dataStr.(field_name).ECL_3sig+0.01]);
xlabel('$$ \bf  Measurement $$', 'Interpreter', 'latex');
ylabel('$$ \bf  ECL\ [eV]  $$', 'Interpreter', 'latex');
ECL_str = string(round(dataStr.(field_name).ECL_mu,3)) + " (" + string(round(dataStr.(field_name).ECL_3sig,3)) + ") eV";
text(0.05,0.8,ECL_str,'Interpreter',"latex",'Units',"normalized",'FontSize',15,'FontWeight','bold');
% Saving the figure
print(path_save + dataStr.(field_name).filename,'-dpng', '-r500');