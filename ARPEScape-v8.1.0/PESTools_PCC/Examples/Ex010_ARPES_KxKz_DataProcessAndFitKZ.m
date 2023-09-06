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
arpes_raw = load_adress_data('010_ARPES_kxkz.h5', path_data);
arpes_dat = delete_scans(arpes_raw, {[2:2:2e3]});   % deleting every even numbered scan, so only the odd-numbered ARPES scans remain
arpes_ref = delete_scans(arpes_raw, {[1:2:2e3]});   % deleting every odd numbered scan, so only the even-numbered XPS scans remain
% -- Viewing the ARPES data
slice_args   = {820, -1.00, 0.00};
view_arpes_data(arpes_dat,slice_args);
view_arpes_data(arpes_ref,slice_args);
%% (1.2) - Binding energy alignment to the reference data
close all; 
% - AUTOMATIC COARSE ALIGNMENT TO THE CORE-LEVEL REFERENCE
alignType   = "align2ref";
scan_indxs	= [];
eWin        = [-2.4, -3];
dEWin       = 1.00;
dESmooth    = [];
feat        = 'peak';
align_args  = {alignType, scan_indxs, eWin, dEWin, dESmooth, feat};
[arpes_dat, arpes_ref]  = align_energy(arpes_dat, align_args, arpes_ref);
view_arpes_data(arpes_dat,slice_args);
view_arpes_data(arpes_ref,slice_args);
% - FINE ALIGNMENT TO KNOW REFERENCE
alignType   = "global shift via scan"; 
scan_indxs	= 30;
eWin        = 2.00;
dEWin       = 1.00;
dESmooth    = [];
feat        = 'edge';
align_args  = {alignType, scan_indxs, eWin, dEWin, dESmooth, feat};
arpes_dat  = align_energy(arpes_dat, align_args);
slice_args   = {820, -0.10, 0.00};
view_arpes_data(arpes_dat,slice_args);
%% (1.3) - Background subtraction and intensity normalisation
close all;
bgrdType    = "AngleInt-clip";
bgrdWin     = [];
bgrdScale   = 0.75;
bgrdSmooth  = 25;
normType    = "max-each";
bgrd_args   = {bgrdType, bgrdWin, bgrdScale, bgrdSmooth, normType}; 
arpes_dat = bgrd_subtr_data(arpes_dat, bgrd_args, 1); 
view_arpes_data(arpes_dat,slice_args);
%% (1.4) - Wave-vector conversions
close all; 
eB_ref      = 0;
kx_ref      = 0;
thtA_ref    = 0;
v000        = 12.57;
inc_alpha   = 9;
kconv_args  = {eB_ref, kx_ref, thtA_ref, v000, inc_alpha}; 
arpes_dat = convert_to_k(arpes_dat, kconv_args, 1);
arpes_dat.kx = arpes_dat.kx + 0.000;
slice_args   = {14.8, -0.10, 0.00};
view_arpes_data(arpes_dat,slice_args);

%% *(2) Luttinger area determination by finding the total area enclosed by the Fermi-surface*
%% (2.1) - Extracting the Fermi-Surface
close all; 
% (A) - Extracting Fermi-Surface Slice
nFrame          = 30;
isoslice_args   = {nFrame, "IsoE", [-0.10, 0.0], 0};
[fs_dat, ~]     = extract_isoSlice(arpes_dat, isoslice_args); fs_dat
% (B) - Background subtracting the Iso-E data
bgrdType    = "AngleInt-Norm";   % either "angle/kx-int (fixed edc)", "angle/kx-int (fixed mean)", "angle/kx-int (moving mean)", "none"
bgrdWin     = [-1.50, -1.00];               % EDC angle range to integrate over to extract ISubtr
bgrdScale   = 0.50;                         % scale coefficient to subtract (Data - intScale * ISubtr)
bgrdSmooth  = 50;                           % single, constant value that determine presmoothing of data to extract background. 
normType    = "none";                  % "init max", "final max", "none"
bgrd_args   = {bgrdType, bgrdWin, bgrdScale, bgrdSmooth, normType}; 
fs_dat    = bgrd_subtr_isoe(fs_dat, bgrd_args); view_arpes_isoSlice(fs_dat);
fs_dat.DSlice = fs_dat.DSlice ./ max(fs_dat.DSlice(:));
%% (2.2) - Extracting iso-contours for number density determination
close all;
icPOI           = [0.00, 14.80];       % 1x2 vector of a given point of interest, aroudn which the contours are determined
icMinArea       = 0.005;            % single, constant value that stores all contours above the minimum enclosed area
icVal           = [0.45, 0.55];     % single value or 1x2 vector that gives the iso-countour value, or range [minVal, maxVal]
icN             = 2;               % single, constant value that determines the total number of iso-contours to extract. For single icVal, this is ignored.
icSmooth        = 2*[1,1];       % Gaussian pre-smoothing parameters of the data prior to extracting the iso-contours
isocont_args    = {nFrame, icPOI, icMinArea, icVal, icN, icSmooth};
[lv_dat, ~] = extract_isoConts(fs_dat, isocont_args); lv_dat
axis([lv_dat.X0-0.25, lv_dat.X0+0.25, lv_dat.Y0-0.25, lv_dat.Y0+0.25]);
lv_dat.Ne
%% (2.3) - Saving corresponding .mat file
close all;
field_name = (filename+"_FitKZ");
% Defining the MATLAB data structure
dataStr.(field_name)                        = struct();
% - Measurement / Sample Information
dataStr.(field_name).filename               = field_name;
dataStr.(field_name).fit_method             = "FitKZ";
dataStr.(field_name).beamtime_id            = beamtime_id;
dataStr.(field_name).dataset_id             = dataset_id;
dataStr.(field_name).sample_id              = sample_id;
dataStr.(field_name).sample_stack_info      = sample_stack_info;
% - Data and Fits
dataStr.(field_name).data                   = arpes_dat;
dataStr.(field_name).fs_dat                 = fs_dat;
dataStr.(field_name).lv_dat                 = lv_dat;
% - Best Fit Results
% -- X0 and Y0 maximum point
dataStr.(field_name).X0                     = lv_dat.X0;
dataStr.(field_name).Y0                     = lv_dat.Y0;
% -- XMag and YMag estimates & measurement uncertainty
dataStr.(field_name).XMag_mu                = lv_dat.XMag(1);
dataStr.(field_name).XMag_3sig              = lv_dat.XMag(2);
dataStr.(field_name).YMag_mu                = lv_dat.YMag(1);
dataStr.(field_name).YMag_3sig              = lv_dat.YMag(2);
% -- NE estimates & measurement uncertainty
dataStr.(field_name).NE_mu                  = lv_dat.Ne(1);
dataStr.(field_name).NE_3sig                = lv_dat.Ne(2);
% - Saving the data
save_adress_data(dataStr.(field_name), path_save + field_name);
%% (2.4) - Saving corresponding .txt file
% Creating an information string
info = "";
info = info + sprintf("Filename: \t\t %s \n",               dataStr.(field_name).filename);
info = info + sprintf("Fit Method: \t %s \n",               dataStr.(field_name).fit_method);
info = info + sprintf("Beamtime ID: \t %s \n",              dataStr.(field_name).beamtime_id);
info = info + sprintf("Dataset ID: \t %s \n",               dataStr.(field_name).dataset_id);
info = info + sprintf("Sample ID: \t\t %s \n",              dataStr.(field_name).sample_id);
info = info + sprintf("Sample stack: \t %s \n",             dataStr.(field_name).sample_stack_info);
info = info + sprintf("NE: \t\t\t %.3s (%.3s) cm^-2 \n",    dataStr.(field_name).NE_mu, dataStr.(field_name).NE_3sig);
info = info + sprintf("X0: \t\t %.3f Ang^-1 \n",            dataStr.(field_name).X0);
info = info + sprintf("Y0: \t\t %.3f Ang^-1 \n",            dataStr.(field_name).Y0);
info = info + sprintf("XMag: \t\t %.3f (%.3f) Ang^-1 \n",   dataStr.(field_name).XMag_mu, dataStr.(field_name).XMag_3sig);
info = info + sprintf("YMag: \t\t %.3f (%.3f) Ang^-1 \n",   dataStr.(field_name).YMag_mu, dataStr.(field_name).YMag_3sig);
info = info + sprintf("Notes: \t\t %s", measurement_notes);
fprintf(info);
% Writing the string to a text file
fileID = fopen(char(path_save + dataStr.(field_name).filename + ".txt"),'w');
fprintf(fileID,info);
fclose(fileID);
%% (2.6) - Saving corresponding .png figure
close all;
pp  = plot_props();
% Initialising the figure
fig = figure();
fig.Position(3) = 3.00*pp.fig4x4(1);
fig.Position(4) = 0.75*pp.fig4x4(2);
% (1) PLOTTING THE KX CUT
scan_index = 30;
kx_slice = extract_isoScan(dataStr.(field_name).data, {scan_index}, 0);
subplot(1,3,1); hold on;
ImData(...
    kx_slice.XScan,...
    kx_slice.YScan,...
    kx_slice.DScan);
img_props(); 
xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
axis([min(kx_slice.XScan(:)), max(kx_slice.XScan(:)), min(kx_slice.YScan(:)), max(kx_slice.YScan(:))]);
title("ARPES: Eb(kx)");
% (2) PLOTTING THE KY CUT
ky_slice = extract_isoSlice(dataStr.(field_name).data, {[], "IsoK", [-0.05, 0.05], 1}, 0);
subplot(1,3,3); hold on;
ImData(...
    ky_slice.XSlice,...
    ky_slice.YSlice,...
    ky_slice.DSlice);
img_props(); 
xlabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
axis([min(ky_slice.XSlice(:)), max(ky_slice.XSlice(:)), min(ky_slice.YSlice(:)), max(ky_slice.YSlice(:))]);
title("ARPES: Eb(kz)");
% (3) PLOTTING THE ISOE SLICE
fs_slice = extract_isoSlice(dataStr.(field_name).data, {[], "IsoE", [-0.10, 0.00], 1}, 0);
subplot(1,3,2); hold on;
% - Plotting the Fermi Surface
ImData(...
    dataStr.(field_name).fs_dat.XSlice,...
    dataStr.(field_name).fs_dat.YSlice,...
    dataStr.(field_name).fs_dat.DSlice);
% - Plotting the iso-contours
iso_cols = jet(length(dataStr.(field_name).lv_dat.XCont));
for i = 1:length(dataStr.(field_name).lv_dat.XCont)
    plot(dataStr.(field_name).lv_dat.XCont{i}, dataStr.(field_name).lv_dat.YCont{i}, 'color', iso_cols(i,:), 'linewidth', 1, 'linestyle', '-');
end
% -- Plotting the total span of the iso-contours
plot(dataStr.(field_name).lv_dat.X0+[-1,1].*0.5.*dataStr.(field_name).lv_dat.XMag(1), dataStr.(field_name).lv_dat.Y0.*[1,1], 'm-', 'linewidth', 1.5);
plot(dataStr.(field_name).lv_dat.X0.*[1,1], dataStr.(field_name).lv_dat.Y0+[-1,1].*0.5.*dataStr.(field_name).lv_dat.YMag(1), 'm-', 'linewidth', 1.5);
axis equal; img_props(); cbar_props(); colorbar off;
xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex');
X0 = dataStr.(field_name).lv_dat.X0;
Y0 = dataStr.(field_name).lv_dat.Y0;
axis([X0-0.35, X0+0.35, Y0-0.35, Y0+0.35]);
xticks(-25:0.1:25); yticks(-25:0.1:25);
title("ARPES: FS(kx,kz)");
print(path_save + dataStr.(field_name).filename,'-dpng', '-r500');
