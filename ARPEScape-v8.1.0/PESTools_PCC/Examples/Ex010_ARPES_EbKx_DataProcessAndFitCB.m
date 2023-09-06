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
end
view_arpes_data(arpes_dat{1});
%% (1.5) - Cropping around the ROI
close all; 
xLims       = [-0.15, 0.15];    % x-axis limits to crop (if empty like [], does not crop this axis)
yLims       = [-0.45, 0.25];     % y-axis limits to crop (if empty like [], does not crop this axis)
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
cutWin          = 0.005*[-1, 1];    % Integration window of the cut to be made.
isocut_args     = {[], cutType, cutWin};
xWin = [];
for i = 1:length(arpes_dat) 
    edc_dat{i}    = extract_isoCut(arpes_dat{i}, isocut_args, 1);
    [edc_dat{i}.peak_loc, edc_dat{i}.peak_int] = find_peak_loc(edc_dat{i}.XCut, edc_dat{i}.DCut, xWin, "spline");
    edc_dat{i}.DCut = edc_dat{i}.DCut ./ edc_dat{i}.peak_int;
end

%% *(2) FitCB2PEAKD2: Find peaks in CB EDC via 2nd derivative*
%% (2.1) - Initialising the ARPES data to be fitted
close all;
arpes_fit               = struct();
arpes_fit.FileName      = file_names;
arpes_fit.repeats       = length(arpes_dat);
arpes_fit.data          = arpes_dat;
arpes_fit.edc           = edc_dat;
%% (2.2) - Extracting peak positions based on 2nd derivative
close all;
EQW = []; arpes_fit.fits = [];
eqwWin = [-0.30, -0.05; -0.05, 0.15];
NQW = size(eqwWin,2);
for i = 1:length(arpes_fit.data)
    fprintf("Run %i / %i \n", i, length(arpes_fit.data));
    % - 1 - Extracting the EDC cut
    XCut    = edc_dat{i}.XCut;
    DCut    = edc_dat{i}.DCut; 
    DCut    = Gaco1(DCut, 1); DCut = DCut - min(DCut(:)); DCut = DCut ./ max(DCut(:)); 
    % - 2 - Finding the second derivative
    d2ydx2  = diff(DCut,2);
    d2ydx2  = Gaco1(d2ydx2, 4); 
    d2ydx2  = d2ydx2 ./ max(d2ydx2(:));
    % - 3 - Finding the second derivative normalised to initial data
    XXCut	= linspace(min(XCut(:)),max(XCut(:)),5e3);
    yy0     = interp1(XCut,DCut,XXCut); 
    yy1     = interp1(XCut(3:end),d2ydx2,XXCut);
    yy1 = yy1 ./ max(yy1(:)); 
    yy2 = yy1 ./ yy0;
    % - 4 - Final data assignment
    yy2 = data_filter1D(yy2, "gaco1",{20});
    DDCut = yy2;
    % - 5 - Cropping around the region of interest
    xWin = [min(eqwWin(:)), max(eqwWin(:))];
    [roi_XCut, roi_DCut]      = data_crop1D(XCut, DCut, xWin);
    [roi_XXCut, roi_DDCut]    = data_crop1D(XXCut, DDCut, xWin);
    % - 6 - Extracting the EQW position from second derivative
    [EQW(:,i), QWS_yval(:,i)] = find_peak_loc(XXCut, -1*DDCut, eqwWin, "spline");
    DDCut = DDCut ./ max(QWS_yval(:,i));
    roi_DDCut = roi_DDCut ./ max(QWS_yval(:,i));
    % - 7 - Final plot of the analysis
    figure(); hold on;
    plot(XCut, DCut, 'k:'); plot(XXCut, DDCut, 'b:');
    plot(roi_XCut, roi_DCut, 'k-', 'linewidth', 1.5); 
    plot(roi_XXCut, roi_DDCut, 'b-', 'linewidth', 1.5);
    for j = 1:NQW
        line([1 1]*EQW(j,i), [-1e5, 1e5], 'Color', [1 0 0], 'LineWidth', 1.5, 'Linestyle', '-');
    end
    gca_props();
    axis([min(XCut(:)), max(XCut(:)), -1.2, 1.2]);
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
    arpes_fit.fits{i}.NQW(:,i)       = NQW;
    arpes_fit.fits{i}.EQW(:,i)       = EQW(:,i);
end
%% (2.3) - Finding the best estimate of EQW across the dataset
close all; 
EQW = [];
for i = 1:length(arpes_fit.fits)
    EQW(:,i) = arpes_fit.fits{i}.EQW(:,i);
end
% -- Isolating any potential outliers
XX              = 1:length(EQW);
TF              = isoutlier(EQW(1,:)); % TF(1:end) = 0;
XX_adj{1}       = XX(~TF);
XX_adj{2}       = XX(TF);
EQW_adj{1}      = EQW(:,~TF);
EQW_adj{2}      = EQW(:,TF);
% -- Plotting the result of outlier analysis
figure(); hold on;
for i = 1:NQW
    plot(XX_adj{1}, EQW_adj{1}(i,:), 'x');
    plot(XX_adj{2}, EQW_adj{2}(i,:), 'ro');
end
% -- Final variable assignments
EQW             = EQW(:,~TF);
EQW_mu          = mean(EQW,2)
EQW_3sig        = 3*std(EQW,[],2)
%% (2.4) - Saving corresponding .mat file
close all;
field_name = (filename+"_FitCB2PEAKD2");
% Defining the MATLAB data structure
dataStr.(field_name)                        = struct();
% - Measurement / Sample Information
dataStr.(field_name).filename               = field_name;
dataStr.(field_name).fit_method             = "FitCB2PEAKD2";
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
% -- EQW estimates & measurement uncertainty
dataStr.(field_name).NQW                    = NQW;
dataStr.(field_name).EQW                    = EQW;
dataStr.(field_name).EQW_mu                 = EQW_mu;
dataStr.(field_name).EQW_3sig               = EQW_3sig;
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
info = info + sprintf("NQW: \t\t\t %.0f \n",  dataStr.(field_name).NQW);
for i = 1:dataStr.(field_name).NQW
    info = info + sprintf("EQW(%.0f): \t\t %.3f (%.3f) eV \n",  i, dataStr.(field_name).EQW_mu(i), dataStr.(field_name).EQW_3sig(i));
end
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
% -- Plotting the best fit and 3sig error of EQW
for i = 1:dataStr.(field_name).NQW
    minX = min(dataStr.(field_name).data{index}.kx(:));
    maxX = max(dataStr.(field_name).data{index}.kx(:));
    patch([minX, maxX, maxX, minX, minX]*1e3,...
        dataStr.(field_name).EQW_mu(i) + [-1, -1, 1, 1, -1]*dataStr.(field_name).EQW_3sig(i),...
        [0 0 0], 'FaceColor', pp.col.fit{i}, 'FaceAlpha', 0.5, 'EdgeAlpha', 0);
    line([minX,maxX]*1e3, [1,1]*dataStr.(field_name).EQW_mu(i), 'Color', pp.col.fit{i}, 'LineWidth', 2, 'Linestyle', '-');
end
% --- Formatting the axis
img_props(); 
cbar_props([], 'Position', [0.06 0.10 0.015 0.175], 'YAxisLocation', 'left');
axis([kx_lims, eb_lims]);
title("ARPES : " + dataStr.(field_name).fit_method);
xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); 
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex'); 

% (2) PLOTTING BEST FIT RESULTS
subplot(2,2,2); hold on;
plot(dataStr.(field_name).fits{index}.XCut, dataStr.(field_name).fits{index}.DCut, 'k:'); 
plot(dataStr.(field_name).fits{index}.XXCut, dataStr.(field_name).fits{index}.DDCut, 'r:');
plot(dataStr.(field_name).fits{index}.roi_XCut, dataStr.(field_name).fits{index}.roi_DCut, 'k-', 'linewidth', 1.5); 
plot(dataStr.(field_name).fits{index}.roi_XXCut, dataStr.(field_name).fits{index}.roi_DDCut, 'r-', 'linewidth', 1.5);
% -- Plotting the best fit and 3sig error of EQW
for i = 1:dataStr.(field_name).NQW
    line([1 1]*EQW_mu(i), [-1e5, 1e5], 'Color', pp.col.fit{i}, 'LineWidth', 2.0, 'Linestyle', '-');
end
gca_props();
axis([min(dataStr.(field_name).fits{index}.XCut(:)), max(dataStr.(field_name).fits{index}.XCut(:)), -1.2, 1.2]);
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Intensity $$', 'Interpreter', 'latex');
% (3) PLOTTING EQW VS MEASUREMENT NUMBER
subplot(2,2,4); hold on;
XX = 1:size(dataStr.(field_name).EQW,2);
for i = 1:dataStr.(field_name).NQW
    plot(XX, dataStr.(field_name).EQW(i,:), 'kx-',...
        'color', pp.col.fit{i}, 'markerfacecolor', pp.col.fit{i}, 'linewidth', 1.5,'markersize', 10);
    errorbar(mean(XX(:)), dataStr.(field_name).EQW_mu(i), dataStr.(field_name).EQW_3sig(i), dataStr.(field_name).EQW_3sig(i), 0, 0, 'ks-',...
        'color', 'k', 'markerfacecolor', pp.col.fit{i}, 'markersize', 10, 'linewidth', 1.5);
    EQW_str = string(round(dataStr.(field_name).EQW_mu(i),3)) + " (" + string(round(dataStr.(field_name).EQW_3sig(i),3)) + ") eV";
    text(1,dataStr.(field_name).EQW_mu(i)+dataStr.(field_name).EQW_3sig(i)+0.005,EQW_str,'Interpreter',"latex",'FontSize',15,'FontWeight','bold');
end
% - Formatting the figure
gca_props(0); 
axis([...
    0.75, max(XX)+0.25,...
    min(dataStr.(field_name).EQW_mu(:)) - 1.25*max(dataStr.(field_name).EQW_3sig(:))-0.01,...
    max(dataStr.(field_name).EQW_mu(:)) + 1.25*max(dataStr.(field_name).EQW_3sig(:))+0.01]);
xlabel('$$ \bf  Measurement $$', 'Interpreter', 'latex');
ylabel('$$ \bf  EQW\ [eV]  $$', 'Interpreter', 'latex');
% Saving the figure
print(path_save + dataStr.(field_name).filename,'-dpng', '-r500');

%% *(3) FitCB2NBOLA: Fit CB data to N-parabolic curves*
%% (3.1) Initialising the ARPES data to be fitted
close all;
arpes_fit               = struct();
arpes_fit.FileName      = file_names;
arpes_fit.repeats       = length(arpes_dat);
arpes_fit.data          = arpes_dat;
%% (3.2) Defining the initial & boundary conditions of best fit
close all;
iparams = {};
% DEFINING THE TYPES OF CURVES TO BE USED
cTYPE   = ["G2DA"; "G2DA"];       % type of curve to use for fitting. Default: "G2DA" ("G2D", "L2D")
% 1 - DEFINING THE INITIAL CONDITIONS OF THE ARPES COMPONENTS
INT     = [0.25; 0.20];          % scalar of the peak intensity of 2D PE curve.
XLOC    = 0.00*[1;1];             % scalar of the x-location of the min/max of the 2D parabolic ARPES dispersion.
YLOC    = [-0.150; -0.030];         % scalar of the y-location of the min/max of the 2D parabolic ARPES dispersion.
XFWHM   = 0.03.*[1;1];         % scalar of the x-axis FWHM for each Gaussian (k-resolution)
YFWHM   = 0.15.*[1;1];         % scalar of the y-axis FWHM for each Gaussian (Eb-resolution)
MSTAR   = 0.06.*[1;1];         % scalar of the effective mass, which governs the curvature of the parabola.
iparams{1} = [INT, XLOC, YLOC, XFWHM, YFWHM, MSTAR]; iparams{1}
% 2 - DEFINING THE UNCERTAINTIES IN THE FIT PARAMETERS
% -- Lower bounds
iparams{2} = iparams{1}; 
iparams{2}(:,1) = iparams{1}(:,1) - 0.10;
iparams{2}(:,2) = iparams{1}(:,2) - 0.002; 
iparams{2}(:,3) = iparams{1}(:,3) - 0.05; 
iparams{2}(:,4) = iparams{1}(:,4) - 0.10; 
iparams{2}(:,5) = iparams{1}(:,5) - 0.10; 
iparams{2}(:,6) = iparams{1}(:,6) - 0.02;
iparams{2}
% -- Upper bounds
iparams{3}      = iparams{1}; 
iparams{3}(:,1) = iparams{1}(:,1) + 0.10;
iparams{3}(:,2) = iparams{1}(:,2) + 0.002; 
iparams{3}(:,3) = iparams{1}(:,3) + 0.05; 
iparams{3}(:,4) = iparams{1}(:,4) + 0.10; 
iparams{3}(:,5) = iparams{1}(:,5) + 0.10; 
iparams{3}(:,6) = iparams{1}(:,6) + 0.02;
iparams{3}
% 3 - DEFINING THE BACKGROUND PARAMETERS
FDEF    =  0.00;        % scalar of the FDD Fermi-Level position.
FDT     =  12.0;        % scalar of the FDD temperature.
FDW     =  0.10;        % scalar of the FDD Gaussian width after convolution.
BGR     =  0.05;        % scalar of the gradient of the linear background.
BIN     =  0.10;        % scalar of the y-intercept of the linear background.
BCO     =  0.05;        % scalar of the constant background y-offset value.
ibgrnd{1} = [FDEF, FDT, FDW, BGR, BIN, BCO]; 
ibgrnd{1}
% 4 - DEFINING THE UNCERTAINTY IN THE BACKGROUND PARAMETERS
% -- Lower bounds
ibgrnd{2}       = ibgrnd{1}; 
ibgrnd{2}(1)    = ibgrnd{1}(1) - 0.00;
ibgrnd{2}(2)    = ibgrnd{1}(2) - 0.00;
ibgrnd{2}(3)    = ibgrnd{1}(3) - 0.10;
ibgrnd{2}(4)    = ibgrnd{1}(4) - 0.10;
ibgrnd{2}(5)    = ibgrnd{1}(5) - 0.10;
ibgrnd{2}(6)    = ibgrnd{1}(6) - 0.10;
ibgrnd{2}
% -- Upper bounds
ibgrnd{3}       = ibgrnd{1}; 
ibgrnd{3}(1)    = ibgrnd{1}(1) + 0.00;
ibgrnd{3}(2)    = ibgrnd{1}(2) + 0.00;
ibgrnd{3}(3)    = ibgrnd{1}(3) + 0.10;
ibgrnd{3}(4)    = ibgrnd{1}(4) + 0.10;
ibgrnd{3}(5)    = ibgrnd{1}(5) + 0.10;
ibgrnd{3}(6)    = ibgrnd{1}(6) + 0.10;
ibgrnd{3}
% 5 - Preview the initial conditions of the model vs data fit
arpes2nbola_view_init(arpes_fit.data{1}, cTYPE, iparams, ibgrnd);
arpes2nbola_view_init(arpes_fit.data{end}, cTYPE, iparams, ibgrnd);
%% (3.3) Execute ARPES fitting algorithm
close all;
solve_type  = "fmincon";
for i = 1:length(arpes_fit.data)
    fprintf("Run %i / %i", i, length(arpes_fit.data));
    arpes_fit.fits{i} = arpes2nbola_solver(arpes_fit.data{i}, cTYPE, iparams, ibgrnd, solve_type);
end
arpes2nbola_view_fit(arpes_fit.fits{1}); arpes_fit.fits{1}
arpes2nbola_view_fit(arpes_fit.fits{end}); arpes_fit.fits{end}
%% (3.4) Finding the best estimate of the band offset
close all;
% 1 - Determine the band offset from the first sub band energy
CHISQ = []; EQW = []; MSTAR = []; 
for i = 1:length(arpes_fit.fits)
    EQW(:,i)    = arpes_fit.fits{i}.YLOC;
    MSTAR(:,i)  = arpes_fit.fits{i}.MSTAR;
    CHISQ(i)    = arpes_fit.fits{i}.CHISQ;
end
% -- Isolating any potential outliers
XX              = 1:length(EQW);
TF              = isoutlier(EQW(1,:)); % TF(1:end) = 0;
XX_adj{1}       = XX(~TF);
XX_adj{2}       = XX(TF);
EQW_adj{1}      = EQW(:,~TF);
EQW_adj{2}      = EQW(:,TF);
% -- Plotting the result of outlier analysis
figure(); hold on;
for i = 1:NQW
    plot(XX_adj{1}, EQW_adj{1}(i,:), 'x');
    plot(XX_adj{2}, EQW_adj{2}(i,:), 'ro');
end
% -- Final variable assignments
EQW            = EQW(:,~TF);
EQW_mu         = mean(EQW,2)
EQW_3sig       = 3*std(EQW,[],2)
CHISQ          = CHISQ(~TF);
CHISQ_mu       = mean(CHISQ(~TF))
CHISQ_3sig     = 3*std(CHISQ(~TF))
MSTAR          = MSTAR(:,~TF);
MSTAR_mu         = mean(MSTAR,2)
MSTAR_3sig       = 3*std(MSTAR,[],2)
%% (3.5) - Saving corresponding .mat file
close all;
field_name = (filename+"_FitCB2NBOLA");
% Defining the MATLAB data structure
dataStr.(field_name)                        = struct();
% - Measurement / Sample Information
dataStr.(field_name).filename               = field_name;
dataStr.(field_name).fit_method             = "FitCB2NBOLA";
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
% -- EQW estimates & measurement uncertainty
dataStr.(field_name).NQW                    = NQW;
dataStr.(field_name).EQW                    = EQW;
dataStr.(field_name).EQW_mu                 = EQW_mu;
dataStr.(field_name).EQW_3sig               = EQW_3sig;
% -- EQW estimates & measurement uncertainty
dataStr.(field_name).MSTAR                  = MSTAR;
dataStr.(field_name).MSTAR_mu               = MSTAR_mu;
dataStr.(field_name).MSTAR_3sig             = MSTAR_3sig;
% -- CHISQ estimates & measurement uncertainty
dataStr.(field_name).CHISQ                  = CHISQ;
dataStr.(field_name).CHISQ_mu               = CHISQ_mu;
dataStr.(field_name).CHISQ_3sig             = CHISQ_3sig;
% - Saving the data
save_adress_data(dataStr.(field_name), path_save + field_name);
%% (3.6) - Saving corresponding .txt file
% Creating an information string
info = "";
info = info + sprintf("Filename: \t\t %s \n",           dataStr.(field_name).filename);
info = info + sprintf("Fit Method: \t %s \n",           dataStr.(field_name).fit_method);
info = info + sprintf("Beamtime ID: \t %s \n",          dataStr.(field_name).beamtime_id);
info = info + sprintf("Dataset ID: \t %s \n",           dataStr.(field_name).dataset_id);
info = info + sprintf("Sample ID: \t\t %s \n",          dataStr.(field_name).sample_id);
info = info + sprintf("Sample stack: \t %s \n",         dataStr.(field_name).sample_stack_info);
info = info + sprintf("Repeats: \t\t %i \n",            dataStr.(field_name).num_of_repeats);
info = info + sprintf("Outliers: \t\t %i \n",           dataStr.(field_name).num_of_outliers);
info = info + sprintf("NQW: \t\t\t %.0f \n",  dataStr.(field_name).NQW);
for i = 1:dataStr.(field_name).NQW
    info = info + sprintf("EQW(%.0f): \t\t %.3f (%.3f) eV \n",  i, dataStr.(field_name).EQW_mu(i), dataStr.(field_name).EQW_3sig(i));
    info = info + sprintf("MSTAR(%.0f): \t\t %.3f (%.3f) \n",  i, dataStr.(field_name).MSTAR_mu(i), dataStr.(field_name).MSTAR_3sig(i));
end
info = info + sprintf("CHISQ: \t\t %.3f (%.3f) \n",  dataStr.(field_name).CHISQ_mu, dataStr.(field_name).CHISQ_3sig);
info = info + sprintf("Notes: \t\t %s", measurement_notes);
fprintf(info);
% Writing the string to a text file
fileID = fopen(char(path_save + dataStr.(field_name).filename + ".txt"),'w');
fprintf(fileID,info);
fclose(fileID);
%% (3.7) - Saving corresponding .png figure
% Initialising the figure
close all;
pp  = plot_props();
fig = figure();
fig.Position(3) = 2.75*pp.fig6x6(1);
fig.Position(4) = 0.65*pp.fig6x6(2);
index = 1;
% (1) PLOTTING THE BEST FITS ON THE ARPES DATA
% -- Plotting the comparitive, final ARPES data
subplot(131); hold on;
% -- Plotting the ARPES data
ImData(...
    dataStr.(field_name).fits{index}.kx,...
    dataStr.(field_name).fits{index}.eb,...
    dataStr.(field_name).fits{index}.D);
% -- Plotting the best fit subbands
for i = 1:dataStr.(field_name).NQW
    plot(dataStr.(field_name).fits{index}.sbn_kx, dataStr.(field_name).fits{index}.sbn_eb{i},...
        '-', 'linewidth', 2, 'color', pp.col.fit{i+1});
end
% --- Formatting the axis
img_props(); 
cbar_props([], 'Position', [0.06 0.10 0.01 0.175], 'YAxisLocation', 'left');
axis([dataStr.(field_name).fits{index}.kx_lims, dataStr.(field_name).fits{index}.eb_lims]);
title('ARPES : ROI Data');
xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); 
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex'); 

% (2) PLOTTING THE BEST FIT TO THE DATA
subplot(132); hold on;
% -- Plotting the ARPES data
ImData(...
    dataStr.(field_name).fits{index}.kx,...
    dataStr.(field_name).fits{index}.eb,...
    dataStr.(field_name).fits{index}.M);
% -- Plotting the best fit subbands
for i = 1:dataStr.(field_name).NQW
    plot(dataStr.(field_name).fits{index}.sbn_kx, dataStr.(field_name).fits{index}.sbn_eb{i},...
        '-', 'linewidth', 2, 'color', pp.col.fit{i+1});
end
% -- Add annotation for the quality of fit
text(0.04, 0.92, "$$ \chi^2 = $$ " + string(dataStr.(field_name).fits{index}.CHISQ),...
    'interpreter', 'latex', 'fontsize', 16, 'color', 'w', 'Units','normalized');
% --- Formatting the axis
img_props();
axis([dataStr.(field_name).fits{index}.kx_lims, dataStr.(field_name).fits{index}.eb_lims]);
title("MODEL : ROI Best Fit : " + dataStr.(field_name).fit_method);
xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex'); 
ylabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex'); 
% (3) PLOTTING EQW VS MEASUREMENT NUMBER
subplot(3,3,3); hold on;
XX = 1:size(dataStr.(field_name).EQW,2);
for i = 1:dataStr.(field_name).NQW
    plot(XX, dataStr.(field_name).EQW(i,:), 'kx-',...
        'color', pp.col.fit{i}, 'markerfacecolor', pp.col.fit{i+1}, 'linewidth', 1.5,'markersize', 10);
    errorbar(mean(XX(:)), dataStr.(field_name).EQW_mu(i), dataStr.(field_name).EQW_3sig(i), dataStr.(field_name).EQW_3sig(i), 0, 0, 'ks-',...
        'color', 'k', 'markerfacecolor', pp.col.fit{i+1}, 'markersize', 10, 'linewidth', 1.5);
    EQW_str = string(round(dataStr.(field_name).EQW_mu(i),3)) + " (" + string(round(dataStr.(field_name).EQW_3sig(i),3)) + ") eV";
    text(1,dataStr.(field_name).EQW_mu(i)+dataStr.(field_name).EQW_3sig(i),EQW_str,'Interpreter',"latex",'FontSize',15,'FontWeight','bold');
end
% - Formatting the figure
gca_props(0); 
axis([...
    0.75, max(XX)+0.25,...
    min(dataStr.(field_name).EQW_mu(:)) - 1.25*max(dataStr.(field_name).EQW_3sig(:))-0.01,...
    max(dataStr.(field_name).EQW_mu(:)) + 1.25*max(dataStr.(field_name).EQW_3sig(:))+0.01]);
ylabel('$$ \bf  EQW\ [eV]  $$', 'Interpreter', 'latex');
% (4) PLOTTING MSTAR VS MEASUREMENT NUMBER
subplot(3,3,6); hold on;
XX = 1:size(dataStr.(field_name).MSTAR,2);
for i = 1:dataStr.(field_name).NQW
    plot(XX, dataStr.(field_name).MSTAR(i,:), 'kx-',...
        'color', pp.col.fit{i}, 'markerfacecolor', pp.col.fit{i+1}, 'linewidth', 1.5,'markersize', 10);
    errorbar(mean(XX(:)), dataStr.(field_name).MSTAR_mu(i), dataStr.(field_name).MSTAR_3sig(i), dataStr.(field_name).MSTAR_3sig(i), 0, 0, 'ks-',...
        'color', 'k', 'markerfacecolor', pp.col.fit{i+1}, 'markersize', 10, 'linewidth', 1.5);
    MSTAR_str = string(round(dataStr.(field_name).MSTAR_mu(i),3)) + " (" + string(round(dataStr.(field_name).MSTAR_3sig(i),3)) + ")";
    text(1,dataStr.(field_name).MSTAR_mu(i)+dataStr.(field_name).MSTAR_3sig(i),MSTAR_str,'Interpreter',"latex",'FontSize',15,'FontWeight','bold');
end
% - Formatting the figure
gca_props(0); 
axis([...
    0.75, max(XX)+0.25,...
    min(dataStr.(field_name).MSTAR_mu(:)) - 1.25*max(dataStr.(field_name).MSTAR_3sig(:))-0.01,...
    max(dataStr.(field_name).MSTAR_mu(:)) + 1.25*max(dataStr.(field_name).MSTAR_3sig(:))+0.01]);
ylabel('$$ \bf  MSTAR\   $$', 'Interpreter', 'latex');
% (5) PLOTTING CHISQ VS MEASUREMENT NUMBER
subplot(3,3,9); hold on;
XX = 1:length(dataStr.(field_name).CHISQ);
plot(XX, dataStr.(field_name).CHISQ, 'kx-',...
    'color', pp.col.fit{1}, 'markerfacecolor', pp.col.fit{1}, 'linewidth', 1.5,'markersize', 10);
errorbar(mean(XX(:)), dataStr.(field_name).CHISQ_mu, dataStr.(field_name).CHISQ_3sig, dataStr.(field_name).CHISQ_3sig, 0, 0, 'ks-',...
    'color', 'k', 'markerfacecolor', pp.col.fit{1}, 'markersize', 10, 'linewidth', 1.5);
CHISQ_str = string(round(dataStr.(field_name).CHISQ_mu,3)) + " (" + string(round(dataStr.(field_name).CHISQ_3sig,3)) + ")";
text(1,dataStr.(field_name).CHISQ_mu+dataStr.(field_name).CHISQ_3sig,CHISQ_str,'Interpreter',"latex",'FontSize',15,'FontWeight','bold');
% - Formatting the figure
gca_props(0); 
axis([...
    0.75, max(XX)+0.25,...
    min(dataStr.(field_name).CHISQ_mu(:)) - 1.25*max(dataStr.(field_name).CHISQ_3sig(:))-0.01,...
    max(dataStr.(field_name).CHISQ_mu(:)) + 1.25*max(dataStr.(field_name).CHISQ_3sig(:))+0.01]);
xlabel('$$ \bf  Measurement $$', 'Interpreter', 'latex');
ylabel('$$ \bf  CHISQ\ [eV]  $$', 'Interpreter', 'latex');
% Saving the figure
print(path_save + dataStr.(field_name).filename,'-dpng', '-r500');