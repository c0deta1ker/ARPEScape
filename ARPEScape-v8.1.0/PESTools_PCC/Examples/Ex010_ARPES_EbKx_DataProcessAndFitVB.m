%% *Processing and plotting ARPES data*
close all; clear all;
%% *Initialising the script*
% Sample & Dataset properties
filename            = "id20230804_ExARPES";
beamtime_id         = "20230804";
dataset_id          = "Example ARPES data";
sample_id           = "('London')";
sample_stack_info   = "";
measurement_notes   = "Sample measured @ADRESS, 77 K.";
% Defining the paths of the model, data and save directories
path_data   = 'D:\OneDrive\PCC_MS_Project\00_SLS_Software\ARPEScape\ARPEScape-v8.0.0\PESTools_PCC\Examples\Example_Data\';
path_save   = 'D:\OneDrive\PCC_MS_Project\00_SLS_Software\ARPEScape\ARPEScape-v8.0.0\PESTools_PCC\Examples\Example_Data_Saved\';

%% *(1) ARPES - Data Processing*
%% (1.1) - Load in data
close all; 
% Load in ARPES data
file_names      = {...
    '010_ARPES_HSCUT.h5',...
    };
arpes_dat = cell(1,length(file_names));
for i = 1:length(file_names); arpes_dat{i} = load_adress_data(file_names{i}, path_data); end
for i = 1:length(arpes_dat); arpes_dat{i} = xcorr_scans(arpes_dat{i}); end
view_arpes_data(arpes_dat{1});
% Load in reference data
file_names_ref      = {...
    '010_ARPES_HSCUT_EF.h5',...
    };
arpes_ref = cell(1,length(file_names_ref));
for i = 1:length(file_names_ref); arpes_ref{i} = load_adress_data(file_names_ref{i}, path_data); end
for i = 1:length(arpes_ref); arpes_ref{i} = xcorr_scans(arpes_ref{i}); end
view_arpes_data(arpes_ref{1});
%% (1.2) - Binding energy alignment to the reference data
close all; 
% - AUTOMATIC COARSE ALIGNMENT TO THE EDGE OF THE REFERENCE
% -- Coarse alignment of the energy
alignType   = "align2ref";
scan_indxs	= [];
eWin        = 0.80;
dEWin       = 1.50;
dESmooth    = [];
feat        = 'edge';
kWin        = 0.95;
align_args  = {alignType, scan_indxs, eWin, dEWin, dESmooth, feat, kWin};
for i = 1:length(arpes_dat)
    [arpes_dat{i}, arpes_ref{i}]  = align_energy(arpes_dat{i}, align_args, arpes_ref{i}, 1);
end
% - FITTING TO THE FDD WITH FINE ALIGNMENT
alignType   = "fit2ef"; 
scan_indxs	= [];
eWin        = 0.00;
dEWin       = 0.25;
dESmooth    = [];
feat        = 'edge';
kWin        = 0.95;
align_args  = {alignType, scan_indxs, eWin, dEWin, dESmooth, feat, kWin};
for i = 1:length(arpes_dat)
    arpes_dat{i}  = align_energy(arpes_dat{i}, align_args, arpes_ref{i}, 1); 
    view_arpes_data(arpes_dat{i}); yline(0, 'color', 'c', 'linewidth', 1.5, 'linestyle', '--'); 
end
%% (1.3) - Background subtraction and intensity normalisation
close all;
bgrdType    = "AngleInt-clip";
bgrdWin     = [];
bgrdScale   = 0.75;
bgrdSmooth  = 25;
normType    = "max-each";
bgrd_args   = {bgrdType, bgrdWin, bgrdScale, bgrdSmooth, normType}; 
for i = 1:length(arpes_dat)
    arpes_dat{i} = bgrd_subtr_data(arpes_dat{i}, bgrd_args, 1); 
end
view_arpes_data(arpes_dat{1});
%% (1.4) - Wave-vector conversions
close all; 
eB_ref      = 0;
kx_ref      = 0;
thtA_ref    = 1.067;
v000        = 12.57;
inc_alpha   = 9;
kconv_args  = {eB_ref, kx_ref, thtA_ref, v000, inc_alpha}; 
for i = 1:length(arpes_dat)
    arpes_dat{i} = convert_to_k(arpes_dat{i}, kconv_args, 1);
    % -- custom shift to centralise kx
    arpes_dat{i}.kx = arpes_dat{i}.kx - 0.0;
end
view_arpes_data(arpes_dat{1});
%% (1.5) - Cropping around the ROI
close all; 
xLims       = [-1.15, 1.15];    % x-axis limits to crop (if empty like [], does not crop this axis)
yLims       = [-6.50, 0.50];     % y-axis limits to crop (if empty like [], does not crop this axis)
zLims       = [];               % z-axis limits to crop (if empty like [], does not crop this axis)
for i = 1:length(arpes_dat) 
    arpes_dat{i}        = data_crop3D(arpes_dat{i}, xLims, yLims, zLims);
    % -- Normalising spectral intensity
    arpes_dat{i}.data   = arpes_dat{i}.data - min(arpes_dat{i}.data(:));
    arpes_dat{i}.data   = arpes_dat{i}.data ./ max(arpes_dat{i}.data(:));
    view_arpes_data(arpes_dat{i}); arpes_dat{i}
end
%% (1.6) - Extracting EDC through k = 0
close all;
edc_dat = {};
cutType         = "EDC";            % Type of cut to make. Either "MDC" (horizontal) or "EDC" (vertical)
cutWin          = 0.01*[-1, 1];    % Integration window of the cut to be made.
isocut_args     = {[], cutType, cutWin};
xWin = [];
for i = 1:length(arpes_dat) 
    edc_dat{i}    = extract_isoCut(arpes_dat{i}, isocut_args, 1);
    [edc_dat{i}.peak_loc, edc_dat{i}.peak_int] = find_peak_loc(edc_dat{i}.XCut, edc_dat{i}.DCut, xWin, "spline");
    edc_dat{i}.DCut = edc_dat{i}.DCut ./ edc_dat{i}.peak_int;
end

%% *(2) FitVB2PEAKD2: Find VBM via 2nd derivative*
%% (2.1) - Initialising the ARPES data to be fitted
close all;
arpes_fit               = struct();
arpes_fit.FileName      = file_names;
arpes_fit.repeats       = length(arpes_dat);
arpes_fit.data          = arpes_dat;
arpes_fit.edc           = edc_dat;
%% (2.2) - Extracting VBM based on 2nd derivative
close all;
VBM = []; arpes_fit.fits = [];
vbmWin = [-2.50, -0.40];
for i = 1:length(arpes_fit.data)
    fprintf("Run %i / %i \n", i, length(arpes_fit.data));
    % - 1 - Extracting the EDC cut
    XCut    = arpes_fit.edc{i}.XCut;
    DCut    = arpes_fit.edc{i}.DCut; 
    DCut    = Gaco1(DCut, 1); DCut = DCut - min(DCut(:)); DCut = DCut ./ max(DCut(:)); 
    % - 2 - Finding the second derivative
    d2ydx2  = diff(DCut,2);
    d2ydx2  = Gaco1(d2ydx2, 20); 
    d2ydx2  = d2ydx2 ./ max(d2ydx2(:));
    % - 3 - Finding the second derivative normalised to initial data
    XXCut	= linspace(min(XCut(:)),max(XCut(:)),5e3);
    yy0     = interp1(XCut,DCut,XXCut); 
    yy1     = interp1(XCut(3:end),d2ydx2,XXCut);
    yy1 = yy1 ./ max(yy1(:)); 
    yy2 = yy1 ./ yy0;
    % - 4 - Final data assignment
    yy2 = data_filter1D(yy2, "gaco1", {20});
    DDCut = yy2;
    % - 5 - Cropping around the region of interest
    xWin = [min(vbmWin(:)), max(vbmWin(:))];
    [roi_XCut, roi_DCut]      = data_crop1D(XCut, DCut, xWin);
    [roi_XXCut, roi_DDCut]    = data_crop1D(XXCut, DDCut, xWin);
    % - 6 - Extracting the VBM position from second derivative
    [VBM, VBM_yval] = find_peak_loc(XXCut, -1*DDCut, vbmWin, "spline");
    DDCut = 0.5*(DDCut ./ max(VBM_yval)+1);
    roi_DDCut = 0.5*(roi_DDCut ./ max(VBM_yval)+1);
    % - 7 - Final plot of the analysis
    figure(); hold on;
    plot(XCut, DCut, 'k-', 'linewidth', 0.5);
    plot(roi_XCut, roi_DCut, 'k-', 'linewidth', 2); 
    plot(XXCut, DDCut, 'b-', 'linewidth', 0.5);
    plot(roi_XXCut, roi_DDCut, 'b-', 'linewidth', 2); 
    a = xline(VBM, 'Color', [1 0 0], 'LineWidth', 1.5, 'Linestyle', '-');
    a.Annotation.LegendInformation.IconDisplayStyle = 'off';
    gca_props();
    axis([min(XCut(:)), max(XCut(:)), min(DCut(:)), max(DCut(:))]);
    gca_props(); title('VBM via 2nd derivative'); 
    xlabel(' X ', 'fontweight', 'bold');
    ylabel(' Y ', 'fontweight', 'bold');
    legend({'Initial Data', 'ROI: Data', 'Initial d2y/dx2', 'ROI: 2y/dx2'}, 'location', 'best', 'fontsize', 9);
    % - 8 - Assigning all variables
    arpes_fit.fits{i}.solve_type = "d2y/dx2 / y";
    % - Assigning the EDCs
    arpes_fit.fits{i}.XCut      = XCut;
    arpes_fit.fits{i}.DCut      = DCut;
    arpes_fit.fits{i}.XXCut     = XXCut;
    arpes_fit.fits{i}.DDCut     = DDCut;
    arpes_fit.fits{i}.roi_XCut  = roi_XCut;
    arpes_fit.fits{i}.roi_DCut  = roi_DCut;
    arpes_fit.fits{i}.roi_XXCut = roi_XXCut;
    arpes_fit.fits{i}.roi_DDCut = roi_DDCut;
    % - Assigning the VBM and BOFF values
    arpes_fit.fits{i}.VBM  = VBM;
end
%% (2.3) - Finding the best estimate of VBM across the dataset
close all; 
VBM = [];
for i = 1:length(arpes_fit.fits)
    VBM(i) = arpes_fit.fits{i}.VBM;
end
% -- Isolating any potential outliers
XX              = 1:length(VBM);
TF              = isoutlier(VBM); % TF(1:end) = 0;
XX_adj{1}       = XX(~TF);
XX_adj{2}       = XX(TF);
VBM_adj{1}      = VBM(:,~TF);
VBM_adj{2}      = VBM(:,TF);
% -- Plotting the result of outlier analysis
figure(); hold on;
plot(XX_adj{1}, VBM_adj{1}, 'kx');
plot(XX_adj{2}, VBM_adj{2}, 'rx');
% -- Final variable assignments
VBM             = VBM(~TF);
VBM_mu          = mean(VBM,2)
VBM_3sig        = 3*std(VBM,[],2)
%% (2.4) - Saving corresponding .mat file
close all;
field_name = (filename+"_FitVB2PEAKD2");
% Defining the MATLAB data structure
dataStr.(field_name)                        = struct();
% - Measurement / Sample Information
dataStr.(field_name).filename               = field_name;
dataStr.(field_name).fit_method             = "FitVB2PEAKD2";
dataStr.(field_name).beamtime_id            = beamtime_id;
dataStr.(field_name).dataset_id             = dataset_id;
dataStr.(field_name).sample_id              = sample_id;
dataStr.(field_name).sample_stack_info      = sample_stack_info;
% - Data and Fits
dataStr.(field_name).num_of_repeats         = length(arpes_fit.fits(~TF));
dataStr.(field_name).num_of_outliers        = length(arpes_fit.fits(TF));
dataStr.(field_name).data                   = arpes_fit.data(~TF);
dataStr.(field_name).edc                    = arpes_fit.edc(~TF);
dataStr.(field_name).fits                   = arpes_fit.fits(~TF);
% - Best Fit Results
% -- VBM estimates & measurement uncertainty
dataStr.(field_name).VBM                    = VBM;
dataStr.(field_name).VBM_mu                 = VBM_mu;
dataStr.(field_name).VBM_3sig               = VBM_3sig;
% - Saving the data
save_adress_data(dataStr.(field_name), path_save + field_name);
%% (2.5) - Saving corresponding .txt file
% Creating an information string
info = "";
info = info + sprintf("Filename: \t\t %s \n",           dataStr.(field_name).filename);
info = info + sprintf("Fit Method: \t %s \n",           dataStr.(field_name).fit_method);
info = info + sprintf("Beamtime ID: \t %s \n",          dataStr.(field_name).beamtime_id);
info = info + sprintf("Dataset ID: \t %s \n",         dataStr.(field_name).dataset_id);
info = info + sprintf("Sample ID: \t\t %s \n",          dataStr.(field_name).sample_id);
info = info + sprintf("Sample stack: \t %s \n",         dataStr.(field_name).sample_stack_info);
info = info + sprintf("Repeats: \t\t %i \n",            dataStr.(field_name).num_of_repeats);
info = info + sprintf("Outliers: \t\t %i \n",           dataStr.(field_name).num_of_outliers);
info = info + sprintf("VBM: \t\t\t %.3f (%.3f) eV \n",  dataStr.(field_name).VBM_mu, dataStr.(field_name).VBM_3sig);
info = info + sprintf("Notes: \t\t %s", measurement_notes);
fprintf(info);
% Writing the string to a text file
fileID = fopen(char(path_save + dataStr.(field_name).filename + ".txt"),'w');
fprintf(fileID,info);
fclose(fileID);
%% (2.6) - Saving corresponding .png figure
close all;
pp  = plot_props();
index = 1;
% Initialising the figure
fig = figure();
fig.Position(3) = 2.00*pp.fig5x5(1);
fig.Position(4) = 0.80*pp.fig5x5(2);

% (1) PLOTTING THE ARPES DATA
kx_lims = [min(dataStr.(field_name).data{index}.kx(:)), max(dataStr.(field_name).data{index}.kx(:))];
eb_lims = [min(dataStr.(field_name).data{index}.eb(:)), max(dataStr.(field_name).data{index}.eb(:))];
% -- Plotting the comparitive, final ARPES data
subplot(1,2,1); hold on;
ImData(...
    dataStr.(field_name).data{index}.kx,...
    dataStr.(field_name).data{index}.eb,...
    dataStr.(field_name).data{index}.data);
% -- Plotting the best fit and 3sig error of VBM
minX = min(dataStr.(field_name).data{index}.kx(:));
maxX = max(dataStr.(field_name).data{index}.kx(:));
patch([minX, maxX, maxX, minX, minX]*1e3,...
    dataStr.(field_name).VBM_mu + [-1, -1, 1, 1, -1]*dataStr.(field_name).VBM_3sig,...
    [0 0 0], 'FaceColor', pp.col.fit{2}, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
yline(dataStr.(field_name).VBM_mu, 'Color', pp.col.fit{2}, 'LineWidth', 2, 'Linestyle', '-');
% --- Formatting the axis
img_props(); 
cbar_props([], 'Position', [0.06 0.10 0.015 0.175], 'YAxisLocation', 'left');
axis([kx_lims, eb_lims]);
title("ARPES : " + dataStr.(field_name).fit_method);
xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); 
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex'); 
% (2) PLOTTING BEST FIT RESULTS
subplot(2,2,2); hold on;
plot(dataStr.(field_name).fits{index}.XCut, dataStr.(field_name).fits{index}.DCut, 'k:', 'linewidth', 0.5); 
plot(dataStr.(field_name).fits{index}.XXCut, dataStr.(field_name).fits{index}.DDCut, 'b:', 'linewidth', 0.5);
plot(dataStr.(field_name).fits{index}.roi_XCut, dataStr.(field_name).fits{index}.roi_DCut, 'k-', 'linewidth', 1); 
plot(dataStr.(field_name).fits{index}.roi_XXCut, dataStr.(field_name).fits{index}.roi_DDCut, 'b-', 'linewidth', 1);
% -- Plotting the best fit and 3sig error of VBM
line([1 1]*VBM_mu, [-1e5, 1e5], 'Color', pp.col.fit{2}, 'LineWidth', 2.0, 'Linestyle', '-');
gca_props();
axis([min(dataStr.(field_name).fits{index}.XCut(:)), max(dataStr.(field_name).fits{index}.XCut(:)), min(dataStr.(field_name).fits{index}.DCut(:)), max(dataStr.(field_name).fits{index}.DCut(:))]);
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Intensity (arb.) $$', 'Interpreter', 'latex');
% (3) PLOTTING VBM VS MEASUREMENT NUMBER
subplot(2,2,4); hold on;
XX = 1:size(dataStr.(field_name).VBM,2);
plot(XX, dataStr.(field_name).VBM, 'kx-',...
    'color', pp.col.fit{i}, 'markerfacecolor', pp.col.fit{2}, 'linewidth', 1.5,'markersize', 10);
errorbar(mean(XX(:)), dataStr.(field_name).VBM_mu, dataStr.(field_name).VBM_3sig, dataStr.(field_name).VBM_3sig, 0, 0, 'ks-',...
    'color', 'k', 'markerfacecolor', pp.col.fit{2}, 'markersize', 10, 'linewidth', 1.5);
VBM_str = string(round(dataStr.(field_name).VBM_mu,3)) + " (" + string(round(dataStr.(field_name).VBM_3sig,3)) + ") eV";
text(1,dataStr.(field_name).VBM_mu+dataStr.(field_name).VBM_3sig+0.005,VBM_str,'Interpreter',"latex",'FontSize',15,'FontWeight','bold');
% - Formatting the figure
gca_props();
axis([...
    0.75, max(XX)+0.25,...
    min(dataStr.(field_name).VBM_mu(:)) - 1.25*max(dataStr.(field_name).VBM_3sig(:))-0.01,...
    max(dataStr.(field_name).VBM_mu(:)) + 1.25*max(dataStr.(field_name).VBM_3sig(:))+0.01]);
xlabel('$$ \bf  Measurement $$', 'Interpreter', 'latex');
ylabel('$$ \bf  VBM\ (eV)  $$', 'Interpreter', 'latex');
% Saving the figure
print(path_save + dataStr.(field_name).filename,'-dpng', '-r500');

%% *(3) FitVB2LE: Find VBM via Leading-Edge Method*
%% (3.1) - Initialising the ARPES data to be fitted
close all;
arpes_fit               = struct();
arpes_fit.FileName      = file_names;
arpes_fit.repeats       = length(arpes_dat);
arpes_fit.data          = arpes_dat;
%% (3.2) - Extracting VBM via Leading-Edge
close all;
arpes_fit.fits       = {};
bgrnd_type      = "poly1";                  % "poly0" (flat background), "poly1" (linear background), "none" (to the x-axis baseline)
kWin_edc        = 0.00 + [-0.02, 0.02];     % window of the EDC taken around the VBM position
eWin_back	    = [0.10, 1.50];          % window of the background region within the EDC cut
eWin_edge	    = [-0.70, -0.60];           % window of the VBM leading edge region within the EDC cut
dESmooth        = [];                       % can be empty []. Single, constant value that defines Gaussian pre-smoothing energy width.
iparams         = {kWin_edc, eWin_back, eWin_edge, dESmooth}; 
for i = 1:length(arpes_dat)
    % -- Finding the VBM position from the leading edge fit
    arpes_fit.fits{i} = vb2le_solver(arpes_dat{i}.kx, arpes_dat{i}.eb, arpes_dat{i}.data, bgrnd_type, iparams, 1)
    title(arpes_dat{i}.FileName, 'Interpreter','none');
end
arpes_fit.fits{1}
%% (3.3) - Finding the best estimate of VBM across the dataset
close all; 
VBM = [];
for i = 1:length(arpes_fit.fits)
    VBM(i) = arpes_fit.fits{i}.VBM;
end
% -- Isolating any potential outliers
XX              = 1:length(VBM);
TF              = isoutlier(VBM); % TF(1:end) = 0;
XX_adj{1}       = XX(~TF);
XX_adj{2}       = XX(TF);
VBM_adj{1}      = VBM(:,~TF);
VBM_adj{2}      = VBM(:,TF);
% -- Plotting the result of outlier analysis
figure(); hold on;
plot(XX_adj{1}, VBM_adj{1}, 'kx');
plot(XX_adj{2}, VBM_adj{2}, 'rx');
% -- Final variable assignments
VBM             = VBM(:,~TF);
VBM_mu          = mean(VBM,2)
VBM_3sig        = 3*std(VBM,[],2)
%% (3.4) - Saving corresponding .mat file
close all;
field_name = (filename+"_FitVB2LE");
% Defining the MATLAB data structure
dataStr.(field_name)                        = struct();
% - Measurement / Sample Information
dataStr.(field_name).filename               = field_name;
dataStr.(field_name).fit_method             = "FitVB2LE";
dataStr.(field_name).beamtime_id            = beamtime_id;
dataStr.(field_name).dataset_id             = dataset_id;
dataStr.(field_name).sample_id              = sample_id;
dataStr.(field_name).sample_stack_info      = sample_stack_info;
% - Data and Fits
dataStr.(field_name).num_of_repeats         = length(arpes_fit.fits(~TF));
dataStr.(field_name).num_of_outliers        = length(arpes_fit.fits(TF));
dataStr.(field_name).data                   = arpes_fit.data(~TF);
dataStr.(field_name).fits                   = arpes_fit.fits(~TF);
% - Best Fit Results
% -- VBM estimates & measurement uncertainty
dataStr.(field_name).VBM                    = VBM;
dataStr.(field_name).VBM_mu                 = VBM_mu;
dataStr.(field_name).VBM_3sig               = VBM_3sig;
% - Saving the data
save_adress_data(dataStr.(field_name), path_save + field_name);
%% (3.5) - Saving corresponding .txt file
% Creating an information string
info = "";
info = info + sprintf("Filename: \t\t %s \n",           dataStr.(field_name).filename);
info = info + sprintf("Fit Method: \t %s \n",           dataStr.(field_name).fit_method);
info = info + sprintf("Beamtime ID: \t %s \n",          dataStr.(field_name).beamtime_id);
info = info + sprintf("Dataset ID: \t %s \n",         dataStr.(field_name).dataset_id);
info = info + sprintf("Sample ID: \t\t %s \n",          dataStr.(field_name).sample_id);
info = info + sprintf("Sample stack: \t %s \n",         dataStr.(field_name).sample_stack_info);
info = info + sprintf("Repeats: \t\t %i \n",            dataStr.(field_name).num_of_repeats);
info = info + sprintf("Outliers: \t\t %i \n",           dataStr.(field_name).num_of_outliers);
info = info + sprintf("VBM: \t\t\t %.3f (%.3f) eV \n",  dataStr.(field_name).VBM_mu, dataStr.(field_name).VBM_3sig);
info = info + sprintf("Notes: \t\t %s", measurement_notes);
fprintf(info);
% Writing the string to a text file
fileID = fopen(char(path_save + dataStr.(field_name).filename + ".txt"),'w');
fprintf(fileID,info);
fclose(fileID);
%% (3.6) - Saving corresponding .png figure
close all;
pp  = plot_props();
index = 1;
% Initialising the figure
fig = figure();
fig.Position(3) = 2.00*pp.fig5x5(1);
fig.Position(4) = 0.80*pp.fig5x5(2);
% (1) PLOTTING THE ARPES DATA
kx_lims = [min(dataStr.(field_name).data{index}.kx(:)), max(dataStr.(field_name).data{index}.kx(:))];
eb_lims = [min(dataStr.(field_name).data{index}.eb(:)), max(dataStr.(field_name).data{index}.eb(:))];
% -- Plotting the comparitive, final ARPES data
subplot(1,2,1); hold on;
ImData(...
    dataStr.(field_name).data{index}.kx,...
    dataStr.(field_name).data{index}.eb,...
    dataStr.(field_name).data{index}.data);
% -- Plotting the best fit and 3sig error of VBM
minX = min(dataStr.(field_name).data{index}.kx(:));
maxX = max(dataStr.(field_name).data{index}.kx(:));
patch([minX, maxX, maxX, minX, minX]*1e3,...
    dataStr.(field_name).VBM_mu + [-1, -1, 1, 1, -1]*dataStr.(field_name).VBM_3sig,...
    [0 0 0], 'FaceColor', pp.col.fit{2}, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
yline(dataStr.(field_name).VBM_mu, 'Color', pp.col.fit{2}, 'LineWidth', 2, 'Linestyle', '-');
% --- Formatting the axis
img_props(); 
cbar_props([], 'Position', [0.06 0.10 0.015 0.175], 'YAxisLocation', 'left');
axis([kx_lims, eb_lims]);
title("ARPES : " + dataStr.(field_name).fit_method);
xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); 
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex'); 
% (2) PLOTTING BEST FIT RESULTS
subplot(2,2,2); hold on;
% -- Plotting the VBM position
line([1 1]*dataStr.(field_name).VBM_mu, [-1e5, 1e5], 'Color', pp.col.fit{2}, 'LineWidth', 2.0, 'Linestyle', '-');
% -- Plotting the full EDC cut first
plot(dataStr.(field_name).fits{index}.XCut, dataStr.(field_name).fits{index}.DCut, 'k.-', 'color', pp.col.fit{1}, 'LineWidth', 0.50);
% -- Plotting the fits to the EDC background
plot(dataStr.(field_name).fits{index}.XCut_back, dataStr.(field_name).fits{index}.DCut_back, 'k.-', 'color', pp.col.fit{2}, 'LineWidth', 1.50);
plot(dataStr.(field_name).fits{index}.XX, dataStr.(field_name).fits{index}.DD_back, 'k:', 'color', pp.col.fit{2}, 'LineWidth', 1.50);
% -- Plotting the fits to the EDC leading edge
plot(dataStr.(field_name).fits{index}.XCut_edge, dataStr.(field_name).fits{index}.DCut_edge, 'g.-', 'color', pp.col.fit{3}, 'LineWidth', 1.50);
plot(dataStr.(field_name).fits{index}.XX, dataStr.(field_name).fits{index}.DD_edge, 'g:', 'color', pp.col.fit{3}, 'LineWidth', 1.50);
% -- Plotting the point of intersection    
plot(dataStr.(field_name).fits{index}.VBM, dataStr.(field_name).fits{index}.VBM_int,...
    'ko', 'markersize', 7, 'color', [0 0 0], 'markerfacecolor', [1 0 0]);
% -- Formatting the figure
gca_props();
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex'); 
ylabel('$$ \bf  Intensity (arb.) $$', 'Interpreter', 'latex');
axis([min(dataStr.(field_name).fits{index}.XCut(:)), max(dataStr.(field_name).fits{index}.XCut(:)),...
    min(dataStr.(field_name).fits{index}.DCut(:)), max(dataStr.(field_name).fits{index}.DCut(:))]);
% (3) PLOTTING VBM VS MEASUREMENT NUMBER
subplot(2,2,4); hold on;
XX = 1:size(dataStr.(field_name).VBM,2);
plot(XX, dataStr.(field_name).VBM, 'kx-',...
    'color', pp.col.fit{1}, 'markerfacecolor', pp.col.fit{2}, 'linewidth', 1.5,'markersize', 10);
errorbar(mean(XX(:)), dataStr.(field_name).VBM_mu, dataStr.(field_name).VBM_3sig, dataStr.(field_name).VBM_3sig, 0, 0, 'ks-',...
    'color', 'k', 'markerfacecolor', pp.col.fit{2}, 'markersize', 10, 'linewidth', 1.5);
VBM_str = string(round(dataStr.(field_name).VBM_mu,3)) + " (" + string(round(dataStr.(field_name).VBM_3sig,3)) + ") eV";
text(1,dataStr.(field_name).VBM_mu+dataStr.(field_name).VBM_3sig+0.005,VBM_str,'Interpreter',"latex",'FontSize',15,'FontWeight','bold');
% - Formatting the figure
gca_props(); 
axis([...
    0.75, max(XX)+0.25,...
    min(dataStr.(field_name).VBM_mu(:)) - 1.25*max(dataStr.(field_name).VBM_3sig(:))-0.01,...
    max(dataStr.(field_name).VBM_mu(:)) + 1.25*max(dataStr.(field_name).VBM_3sig(:))+0.01]);
xlabel('$$ \bf  Measurement $$', 'Interpreter', 'latex');
ylabel('$$ \bf  VBM\ (eV)  $$', 'Interpreter', 'latex');
% Saving the figure
print(path_save + dataStr.(field_name).filename,'-dpng', '-r500');

%% *(4) FitVB2LP: Find VBM via Leading-Peak Method*
%% (4.1) - Initialising the ARPES data to be fitted
close all;
arpes_fit               = struct();
arpes_fit.FileName      = file_names;
arpes_fit.repeats       = length(arpes_dat);
arpes_fit.data          = arpes_dat;
%% (4.2) - Extracting VBM via Leading-Edge
close all;
cutType         = "EDC";                    % Type of cut to make. Either "MDC" (horizontal) or "EDC" (vertical)
cutWin          = 0.00 + [-0.02, 0.02];     % Integration window of the cut to be made.
isocut_args     = {[], cutType, cutWin};
xWin         = [-0.90, -0.55];
for i = 1:length(arpes_fit.data) 
    arpes_fit.cuts{i} = extract_isoCut(arpes_dat{i}, isocut_args, 1);
    [arpes_fit.VBM(i), ~] = find_peak_loc(arpes_fit.cuts{i}.XCut, arpes_fit.cuts{i}.DCut, xWin, "sGLA", 1);
end
%% (4.3) - Finding the best estimate of VBM across the dataset
VBM     = arpes_fit.VBM;
% -- Isolating any potential outliers
XX              = 1:length(VBM);
TF              = isoutlier(VBM); TF(1:end) = 0;
XX_adj{1}       = XX(~TF);
XX_adj{2}       = XX(TF);
VBM_adj{1}      = VBM(~TF);
VBM_adj{2}      = VBM(TF);
% -- Redefining arpes_fit.fits and VBM
arpes_fit.cuts  = arpes_fit.cuts(~TF);
VBM             = VBM(~TF);
VBM_mu          = mean(VBM)
VBM_3sig        = 3*std(VBM)
%% (4.4) - Saving corresponding .mat file
close all;
field_name = (filename+"_FitVB2LP");
% Defining the MATLAB data structure
dataStr.(field_name)                        = struct();
% - Measurement / Sample Information
dataStr.(field_name).filename               = field_name;
dataStr.(field_name).fit_method             = "FitVB2LP";
dataStr.(field_name).beamtime_id            = beamtime_id;
dataStr.(field_name).dataset_id             = dataset_id;
dataStr.(field_name).sample_id              = sample_id;
dataStr.(field_name).sample_stack_info      = sample_stack_info;
% - Data and Fits
dataStr.(field_name).num_of_repeats         = length(arpes_fit.cuts(~TF));
dataStr.(field_name).num_of_outliers        = length(arpes_fit.cuts(TF));
dataStr.(field_name).data                   = arpes_fit.data(~TF);
dataStr.(field_name).cuts                   = arpes_fit.cuts(~TF);
% - Best Fit Results
% -- VBM estimates & measurement uncertainty
dataStr.(field_name).VBM                    = VBM;
dataStr.(field_name).VBM_mu                 = VBM_mu;
dataStr.(field_name).VBM_3sig               = VBM_3sig;
% - Saving the data
save_adress_data(dataStr.(field_name), path_save + field_name);
%% (4.5) - Saving corresponding .txt file
% Creating an information string
info = "";
info = info + sprintf("Filename: \t\t %s \n",           dataStr.(field_name).filename);
info = info + sprintf("Fit Method: \t %s \n",           dataStr.(field_name).fit_method);
info = info + sprintf("Beamtime ID: \t %s \n",          dataStr.(field_name).beamtime_id);
info = info + sprintf("Dataset ID: \t %s \n",         dataStr.(field_name).dataset_id);
info = info + sprintf("Sample ID: \t\t %s \n",          dataStr.(field_name).sample_id);
info = info + sprintf("Sample stack: \t %s \n",         dataStr.(field_name).sample_stack_info);
info = info + sprintf("Repeats: \t\t %i \n",            dataStr.(field_name).num_of_repeats);
info = info + sprintf("Outliers: \t\t %i \n",           dataStr.(field_name).num_of_outliers);
info = info + sprintf("VBM: \t\t\t %.3f (%.3f) eV \n",  dataStr.(field_name).VBM_mu, dataStr.(field_name).VBM_3sig);
info = info + sprintf("Notes: \t\t %s", measurement_notes);
fprintf(info);
% Writing the string to a text file
fileID = fopen(char(path_save + dataStr.(field_name).filename + ".txt"),'w');
fprintf(fileID,info);
fclose(fileID);
%% (4.6) - Saving corresponding .png figure
close all;
pp  = plot_props();
index = 1;
% Initialising the figure
fig = figure();
fig.Position(3) = 2.00*pp.fig5x5(1);
fig.Position(4) = 0.80*pp.fig5x5(2);
% (1) PLOTTING THE ARPES DATA
kx_lims = [min(dataStr.(field_name).data{index}.kx(:)), max(dataStr.(field_name).data{index}.kx(:))];
eb_lims = [min(dataStr.(field_name).data{index}.eb(:)), max(dataStr.(field_name).data{index}.eb(:))];
% -- Plotting the comparitive, final ARPES data
subplot(1,2,1); hold on;
ImData(...
    dataStr.(field_name).data{index}.kx,...
    dataStr.(field_name).data{index}.eb,...
    dataStr.(field_name).data{index}.data);
% -- Plotting the best fit and 3sig error of VBM
minX = min(dataStr.(field_name).data{index}.kx(:));
maxX = max(dataStr.(field_name).data{index}.kx(:));
patch([minX, maxX, maxX, minX, minX]*1e3,...
    dataStr.(field_name).VBM_mu + [-1, -1, 1, 1, -1]*dataStr.(field_name).VBM_3sig,...
    [0 0 0], 'FaceColor', pp.col.fit{2}, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
yline(dataStr.(field_name).VBM_mu, 'Color', pp.col.fit{2}, 'LineWidth', 2, 'Linestyle', '-');
% --- Formatting the axis
img_props(); 
cbar_props([], 'Position', [0.06 0.10 0.015 0.175], 'YAxisLocation', 'left');
axis([kx_lims, eb_lims]);
title("ARPES : " + dataStr.(field_name).fit_method);
xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); 
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex'); 

% (2) PLOTTING BEST FIT RESULTS
subplot(2,2,2); hold on;
% -- Plotting the VBM position
line([1 1]*dataStr.(field_name).VBM_mu, [-1e5, 1e5], 'Color', pp.col.fit{2}, 'LineWidth', 2.0, 'Linestyle', '-');
% -- Plotting the full EDC cut first
plot(dataStr.(field_name).cuts{index}.XCut, dataStr.(field_name).cuts{index}.DCut, 'k.-', 'color', pp.col.fit{1}, 'LineWidth', 0.50);
% -- Formatting the figure
gca_props();
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex'); 
ylabel('$$ \bf  Intensity (arb.) $$', 'Interpreter', 'latex');
axis([min(dataStr.(field_name).cuts{index}.XCut(:)), max(dataStr.(field_name).cuts{index}.XCut(:)),...
    min(dataStr.(field_name).cuts{index}.DCut(:)), max(dataStr.(field_name).cuts{index}.DCut(:))]);

% (3) PLOTTING VBM VS MEASUREMENT NUMBER
subplot(2,2,4); hold on;
XX = 1:size(dataStr.(field_name).VBM,2);
plot(XX, dataStr.(field_name).VBM, 'kx-',...
    'color', pp.col.fit{1}, 'markerfacecolor', pp.col.fit{2}, 'linewidth', 1.5,'markersize', 10);
errorbar(mean(XX(:)), dataStr.(field_name).VBM_mu, dataStr.(field_name).VBM_3sig, dataStr.(field_name).VBM_3sig, 0, 0, 'ks-',...
    'color', 'k', 'markerfacecolor', pp.col.fit{2}, 'markersize', 10, 'linewidth', 1.5);
VBM_str = string(round(dataStr.(field_name).VBM_mu,3)) + " (" + string(round(dataStr.(field_name).VBM_3sig,3)) + ") eV";
text(1,dataStr.(field_name).VBM_mu+dataStr.(field_name).VBM_3sig+0.005,VBM_str,'Interpreter',"latex",'FontSize',15,'FontWeight','bold');
% - Formatting the figure
gca_props(); 
axis([...
    0.75, max(XX)+0.25,...
    min(dataStr.(field_name).VBM_mu(:)) - 1.25*max(dataStr.(field_name).VBM_3sig(:))-0.01,...
    max(dataStr.(field_name).VBM_mu(:)) + 1.25*max(dataStr.(field_name).VBM_3sig(:))+0.01]);
xlabel('$$ \bf  Measurement $$', 'Interpreter', 'latex');
ylabel('$$ \bf  VBM\ (eV)  $$', 'Interpreter', 'latex');
% Saving the figure
print(path_save + dataStr.(field_name).filename,'-dpng', '-r500');