%% *Initialising the script*
close all; clear all;
% Sample & Dataset properties
filename            = "id20230806_ExXPS";
beamtime_id         = "20230806";
dataset_id          = "Example ARPES data";
sample_id           = "('London')";
sample_stack_info   = "";
measurement_notes   = "Sample measured @ADRESS, 77 K.";
% Defining the paths of the model, data and save directories
path_data   = 'D:\OneDrive\PCC_MS_Project\00_SLS_Software\ARPEScape\ARPEScape-v8.0.0\PESTools_PCC\Examples\Example_Data\';
path_save   = 'D:\OneDrive\PCC_MS_Project\00_SLS_Software\ARPEScape\ARPEScape-v8.0.0\PESTools_PCC\Examples\Example_Data_Saved\';

%% *(1) XPS - Data Processing*
%% (1.1) - Load in data
close all; 
% Load in XPS data
file_names      = {...
    '020_XPS_In4d_hv=550eV.h5',...
    '020_XPS_In4d_hv=750eV.h5',...
    '020_XPS_In4d_hv=950eV.h5',...
    '020_XPS_In4d_hv=1150eV.h5',...
    };
xps_dat = cell(1,length(file_names));
for i = 1:length(file_names); xps_dat{i} = load_adress_data(file_names{i}, path_data); end
view_arpes_data(xps_dat{1});
% Load in Reference data
file_names_ef      = {...
    '020_XPS_In4d_hv=550eV_EF.h5',...
    '020_XPS_In4d_hv=750eV_EF.h5',...
    '020_XPS_In4d_hv=950eV_EF.h5',...
    '020_XPS_In4d_hv=1150eV_EF.h5',...
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
eWin        = 0.00;
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

%% *(2) XPS - In4d versus photon energy*
%% (2.1) - XPS Curve Fitting (Initialisation)
close all;
% - Initialising data structure
xps_fit               = struct();
xps_fit.FileName      = filename;
xps_fit.repeats       = length(xps_dat);
xps_fit.data          = xps_dat;
for i = 1:length(xps_fit.data); xps_fit.hv(i) = round(xps_fit.data{i}.hv,0); end
% DEFINING THE TYPES OF CURVES TO BE USED
cTYPE   = ["sGLA";  "sGLA";  "sGLA";  "sGLA"];  % type of curve to use for fitting. Default: "sGLA" ("pGLA", "DS")
% 1.1 - DEFINING THE INITIAL CONDITIONS OF THE XPS COMPONENTS
BE      = [-17.47;  -17.95;  -16.77; -19.32]; % scalar of the binding energy of PE curve. Each new row gives a new BE component. Make sure sizes are all consistent.
INT     = [0.70;     0.08;   0.18;  0.34];    % scalar of the peak intensity of PE curve.
FWHM    = [0.38;    0.45;   0.49;   0.31];    % scalar of the FWHM of the PE curve.
MR      = [0.50;    0.50;   0.50;   0.50];    % scalar of the Mixing Ratio: 0 for pure Gaussian, 1 is for pure Lorentzian.
LSE     = [-0.85;  -0.85;   -0.85;  -0.42];   % scalar of the binding energy of spin-orbit split PE curve.
LSI     = [0.66;    0.66;   0.66;   0.55];    % scalar of the branching ratio of spin-orbit split PE curve.
LSW     = [0.0;     0.0;    0.0;    0.0];     % scalar of the additional lorentzian width of spin-orbit split PE curve.
ASY     = [0.0;     0.0;    0.0;    0.0];     % scalar of the PE curve asymmetry parameter (usually for metallic systems).
iparams{1} = [BE, INT, FWHM, MR, LSE, LSI, LSW, ASY]; 
iparams{1}
% 1.2 - DEFINING THE UNCERTAINTIES IN THE FIT PARAMETERS
% -- Lower bounds
iparams{2} = iparams{1}; 
iparams{2}(:,1) = iparams{1}(:,1) - 0.20;
iparams{2}(:,2) = iparams{1}(:,2) - 0.75;
iparams{2}(:,3) = iparams{1}(:,3) - 0.30; 
iparams{2}(:,4) = iparams{1}(:,4) - 0.25; 
iparams{2}(:,5) = iparams{1}(:,5) - 0.05; 
iparams{2}(:,6) = iparams{1}(:,6) - 0.10;  
iparams{2}(:,7) = iparams{1}(:,7) - 0.00; 
iparams{2}(:,8) = iparams{1}(:,8) - 0.00;
iparams{2}
% -- Upper bounds
iparams{3}      = iparams{1}; 
iparams{3}(:,1) = iparams{1}(:,1) + 0.20; 
iparams{3}(:,2) = iparams{1}(:,2) + 0.75;
iparams{3}(:,3) = iparams{1}(:,3) + 0.30; 
iparams{3}(:,4) = iparams{1}(:,4) + 0.25; 
iparams{3}(:,5) = iparams{1}(:,5) + 0.05; 
iparams{3}(:,6) = iparams{1}(:,6) + 0.10; 
iparams{3}(:,7) = iparams{1}(:,7) + 0.00; 
iparams{3}(:,8) = iparams{1}(:,8) + 0.00;
iparams{3}
% -- Forced constraints
iparams{2}(3,6) = iparams{1}(3,6) - 0.0; iparams{3}(3,6) = iparams{1}(3,6) - 0.0;

% 2.1 - DEFINING THE BACKGROUND PARAMETERS
help PESBackground;
bTYPE   = "poly";       % type of background to use for fitting. Default: "Poly" ("none", "Shir", "LinShir")
LHS     = -20.5;      % scalar of the START point on the LHS for background
RHS     = -15.5;      % scalar of the END point on the RHS for background
BGR     = 0.00;      % scalar for a constant background to be included in the fit
ibgrnd{1} = [LHS, RHS, BGR]; 
ibgrnd{1}
% 2.2 - DEFINING THE UNCERTAINTY IN THE BACKGROUND PARAMETERS
% -- Lower bounds
ibgrnd{2} = abs(0.0.*ibgrnd{1}); 
ibgrnd{2}(:,1) = ibgrnd{1}(:,1) - 0.00;
ibgrnd{2}(:,2) = ibgrnd{1}(:,2) - 0.00;
ibgrnd{2}(:,3) = ibgrnd{1}(:,3) - 0.10;
ibgrnd{2}
% -- Upper bounds
ibgrnd{3} = abs(0.0.*ibgrnd{1}); 
ibgrnd{3}(:,1) = ibgrnd{1}(:,1) + 0.00;
ibgrnd{3}(:,2) = ibgrnd{1}(:,2) + 0.00;
ibgrnd{3}(:,3) = ibgrnd{1}(:,3) + 0.10;
ibgrnd{3}
% 2.3 - DEFINING THE BACKGROUND ARGUMENTS
ibgrnd{4} = {[2]};

% 3 - Preview the initial conditions of the model vs data fit
pes2ncurve_view_init(xps_fit.data{1}, cTYPE, iparams, bTYPE, ibgrnd);
pes2ncurve_view_init(xps_fit.data{end}, cTYPE, iparams, bTYPE, ibgrnd);

%% (2.3) - XPS Curve Fitting (Execution)
close all;
solve_type  = "fmincon";
for i = 1:length(xps_fit.data)
    fprintf("Run %i / %i", i, length(xps_fit.data));
    % -- Solve the fit
    xps_fit.fits{i} = pes2ncurve_solver(xps_fit.data{i}, cTYPE, iparams, bTYPE, ibgrnd, solve_type);
    % -- Saving the figure
    pes2ncurve_view_fit(xps_fit.fits{i}); print(path_save + filename + "_" + string(xps_fit.data{i}.hv) + "eV",'-dpng', '-r500');
end
pes2ncurve_view_fit(xps_fit.fits{1});
pes2ncurve_view_fit(xps_fit.fits{end});
%% (2.4) - Extracting best fit parameters for all fits
BE = []; AREA0 = []; CHISQ = [];
for i = 1:length(xps_fit.fits)
    BE(:,i)         = xps_fit.fits{i}.BE;
    AREA0(:,i)      = xps_fit.fits{i}.AREA0;
    CHISQ(:,i)      = xps_fit.fits{i}.CHISQ;
end
% -- Extracting mean and 3-sig values
BE_mu       = mean(BE,2)
BE_3sig     = 3*std(BE,0,2)
AREA0_mu    = mean(AREA0,2)
AREA0_3sig  = 3*std(AREA0,0,2)
CHISQ_mu    = mean(CHISQ,2)
CHISQ_3sig  = 3*std(CHISQ,0,2)
%% (2.5) - Saving corresponding .mat file
field_name = (filename);
% Defining the MATLAB data structure
dataStr.(field_name)                        = struct();
% - Measurement / Sample Information
dataStr.(field_name).filename               = field_name;
dataStr.(field_name).fit_method             = "DepthProfileXPS";
dataStr.(field_name).beamtime_id            = beamtime_id;
dataStr.(field_name).dataset_id             = dataset_id;
dataStr.(field_name).sample_id              = sample_id;
dataStr.(field_name).sample_stack_info      = sample_stack_info;
dataStr.(field_name).measurement_notes      = measurement_notes;
% - Best Fit XPS Results
% -- XPS data and corresponding best fits
dataStr.(field_name).num_of_repeats         = length(xps_fit.data);
dataStr.(field_name).xps_hv                 = xps_fit.hv;
dataStr.(field_name).xps_data               = xps_fit.data;
dataStr.(field_name).xps_fits               = xps_fit.fits;
% -- BE estimates & measurement uncertainty
dataStr.(field_name).BE                    = BE;
dataStr.(field_name).BE_mu                 = BE_mu;
dataStr.(field_name).BE_3sig               = BE_3sig;
% -- AREA0 estimates & measurement uncertainty
dataStr.(field_name).AREA0                 = AREA0;
dataStr.(field_name).AREA0_mu              = AREA0_mu;
dataStr.(field_name).AREA0_3sig            = AREA0_3sig;
% -- CHISQ estimates & measurement uncertainty
dataStr.(field_name).CHISQ                 = CHISQ;
dataStr.(field_name).CHISQ_mu              = CHISQ_mu;
dataStr.(field_name).CHISQ_3sig            = CHISQ_3sig;
% - Saving the data
save_adress_data(dataStr.(field_name), path_save + field_name);
%% (2.6) - Saving corresponding .txt file
% Creating an information string
info = "";
info = info + sprintf("Filename: \t\t %s \n",   dataStr.(field_name).filename);
info = info + sprintf("Fit Method: \t %s \n",   dataStr.(field_name).fit_method);
info = info + sprintf("Beamtime ID: \t %s \n",  dataStr.(field_name).beamtime_id);
info = info + sprintf("Dataset ID: \t\t %s \n",  dataStr.(field_name).dataset_id);
info = info + sprintf("Sample ID: \t\t %s \n",  dataStr.(field_name).sample_id);
info = info + sprintf("Sample stack: \t %s \n", dataStr.(field_name).sample_stack_info);
info = info + sprintf("Repeats: \t\t %i \n",    dataStr.(field_name).num_of_repeats);
info = info + sprintf("BE(InAs): \t\t\t %.3f (%.3f) eV \n", dataStr.(field_name).BE_mu(1), dataStr.(field_name).BE_3sig(1));
info = info + sprintf("BE(InAs(SCLS)): \t\t\t %.3f (%.3f) eV \n", dataStr.(field_name).BE_mu(2), dataStr.(field_name).BE_3sig(2));
info = info + sprintf("CHISQ: \t\t\t %.3f (%.3f) \n", dataStr.(field_name).CHISQ_mu, dataStr.(field_name).CHISQ_3sig);
info = info + sprintf("Notes: \t\t %s", measurement_notes);
fprintf(info);
% Writing the string to a text file
fileID = fopen(char(path_save + dataStr.(field_name).filename + ".txt"),'w');
fprintf(fileID,info);
fclose(fileID);
%% (2.7) - Saving corresponding figure summary
close all;
pp  = plot_props();
index = 1;
% Initialising the figure
fig = figure();
fig.Position(3) = 3.50*pp.fig6x6(1);
fig.Position(4) = 1.25*pp.fig6x6(2);

% Plotting all of the XPS data
subplot(4,3,[1,4]); hold on;
cols = parula(length(dataStr.(field_name).xps_data));
x_vals = []; y_vals = [];
for i = 1:length(dataStr.(field_name).xps_data)
    plot(dataStr.(field_name).xps_data{i}.xdat, dataStr.(field_name).xps_data{i}.ydat,...
        'k-', 'color', cols(i,:), 'linewidth', 2);
    x_vals = vertcat(x_vals, dataStr.(field_name).xps_data{i}.xdat);
    y_vals = vertcat(y_vals, dataStr.(field_name).xps_data{i}.ydat);
end
% - Formatting the figure
gca_props(); % axis tight;
% xlabel('$$ \bf E_B\ -\ E_F\ (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf Intensity\ (arb.) $$', 'Interpreter', 'latex');
axis([min(x_vals(:)), max(x_vals(:)), min(y_vals(:)), 1.1*max(y_vals(:))]);
title_text = dataStr.(field_name).sample_id + " : All Data";
title(title_text, 'interpreter', 'none');

% Plotting the XPS data and background subtraction
subplot(4,3,[7,10]); hold on;
% - Plotting the ROI, LHS and RHS analysis windows
ROI     = [dataStr.(field_name).xps_fits{index}.X(1), dataStr.(field_name).xps_fits{index}.X(end)]; 
LHS     = ROI(1) + 0.05.*[-1,1].*range(dataStr.(field_name).xps_fits{index}.xdat(:));
RHS     = ROI(2) + 0.05.*[-1,1].*range(dataStr.(field_name).xps_fits{index}.xdat(:));
hLHS    = patch([LHS(1), LHS(1), LHS(2), LHS(2), LHS(1)], [-1, 1, 1, -1, -1].*1e6, [0.6 0.6 0.2], 'facealpha', 0.4, 'edgecolor', 'none');
hLHS.Annotation.LegendInformation.IconDisplayStyle = 'off';
hRHS    = patch([RHS(1), RHS(1), RHS(2), RHS(2), RHS(1)], [-1, 1, 1, -1, -1].*1e6, [0.6 0.6 0.2], 'facealpha', 0.4, 'edgecolor', 'none');
hRHS.Annotation.LegendInformation.IconDisplayStyle = 'off';
hMID    = patch([ROI(1), ROI(1), ROI(2), ROI(2), ROI(1)], [-1, 1, 1, -1, -1].*1e6, [0.8 0.9 0.8], 'facealpha', 0.5, 'edgecolor', 'none');
hMID.Annotation.LegendInformation.IconDisplayStyle = 'off';
% - Plotting a vertical line to show the ROI
a = line([ROI(1) ROI(1)], [-1e5, 1e5], 'Color', [0 0 0], 'LineWidth', 1, 'Linestyle', '-');
a.Annotation.LegendInformation.IconDisplayStyle = 'off';
b = line([ROI(2) ROI(2)], [-1e5, 1e5], 'Color', [0 0 0], 'LineWidth', 1, 'Linestyle', '-');
b.Annotation.LegendInformation.IconDisplayStyle = 'off';
% - Plotting the 1D data
plot(dataStr.(field_name).xps_fits{index}.xdat, dataStr.(field_name).xps_fits{index}.ydat, 'b-', 'linewidth', 0.5);
plot(dataStr.(field_name).xps_fits{index}.X, dataStr.(field_name).xps_fits{index}.D, 'b-', 'linewidth', 2);
plot(dataStr.(field_name).xps_fits{index}.X, dataStr.(field_name).xps_fits{index}.B, 'r-', 'linewidth', 2);
plot(dataStr.(field_name).xps_fits{index}.X, dataStr.(field_name).xps_fits{index}.DB, 'k-', 'linewidth', 2);
gca_props(); 
xlabel('$$ \bf  E_B - E_F (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Intensity$$', 'Interpreter', 'latex');
title('Background Subtraction', 'interpreter', 'none'); 
legend({'Initial Data', 'ROI: Data', 'ROI: Background', 'ROI: Final'}, 'location', 'best', 'fontsize', 9);
% - Determining the best limits for the plot
axLim_y = [dataStr.(field_name).xps_fits{index}.ydat; dataStr.(field_name).xps_fits{index}.DB; dataStr.(field_name).xps_fits{index}.B; dataStr.(field_name).xps_fits{index}.D];
axis([...
    min(dataStr.(field_name).xps_fits{index}.xdat(:)), max(dataStr.(field_name).xps_fits{index}.xdat(:)),...
    min(axLim_y(:)), 1.1*max(axLim_y(:))]);

% Plotting best fit XPS data
subplot(4,3,[2,5,8]); hold on;
% -- Plotting all of the curve components
for i = 1:dataStr.(field_name).xps_fits{index}.nSTATES
    area(dataStr.(field_name).xps_fits{index}.XX, dataStr.(field_name).xps_fits{index}.cYY(:,i), 'FaceColor', pp.col.fit{i}, 'FaceAlpha', pp.falpha, 'EdgeAlpha', 0);
    plot(dataStr.(field_name).xps_fits{index}.XX, dataStr.(field_name).xps_fits{index}.cYY(:,i), 'k-', 'linewidth', 0.25);
end
% -- Plotting the curve component energy locations
for i = 1:dataStr.(field_name).xps_fits{index}.nSTATES
    if i == 1; lWidth = pp.llwidth; else; lWidth = pp.lwidth; end
    % --- Plotting BE
    line([1,1]*dataStr.(field_name).xps_fits{index}.BE(i), [0, max(dataStr.(field_name).xps_fits{index}.cYY(:,i))], 'Color', [0 0 0], 'LineWidth', lWidth, 'Linestyle', '-');
end
% -- Plotting the experimental and fit spectra
plot(dataStr.(field_name).xps_fits{index}.X, dataStr.(field_name).xps_fits{index}.DB, 'k-', 'color', pp.col.dat{1}, 'linewidth', 2*pp.llwidth);
plot(dataStr.(field_name).xps_fits{index}.X, dataStr.(field_name).xps_fits{index}.M, 'k-', 'Color', pp.col.dat{2}, 'linewidth', 2*pp.lwidth);
% -- Add annotation for the quality of fit
text(0.04, 0.90, "$$ \chi^2 = $$ " + string(dataStr.(field_name).xps_fits{index}.CHISQ),...
    'interpreter', 'latex', 'fontsize', 15, 'color', 'k', 'Units','normalized');
% - Formatting the figure
gca_props(); % grid on;
ylabel('$$ \bf  Intensity$$', 'Interpreter', 'latex');
axLim_y = [0; dataStr.(field_name).xps_fits{index}.M; dataStr.(field_name).xps_fits{index}.DB; dataStr.(field_name).xps_fits{index}.cYY(:)];
axis([...
    min(dataStr.(field_name).xps_fits{index}.X(:)), max(dataStr.(field_name).xps_fits{index}.X(:)),...
    min(axLim_y(:)), 1.10*max(axLim_y(:))]);
title("Best Curve Fit", 'interpreter', 'none');
% Plotting the residuals
subplot(4,3,[11]); hold on;
bar(dataStr.(field_name).xps_fits{index}.X, dataStr.(field_name).xps_fits{index}.R, 'facecolor', [0 0 0]);
gca_props(); % grid on;
xlabel('$$ \bf  E_B\ -\ E_F\ (eV) $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Resid. $$', 'Interpreter', 'latex');
axLim_y = abs(dataStr.(field_name).xps_fits{index}.R);
axis([...
    min(dataStr.(field_name).xps_fits{index}.X(:)), max(dataStr.(field_name).xps_fits{index}.X(:)),...
    -1.10*max(axLim_y(:)), 1.10*max(axLim_y(:))]);

% Plotting AREA0
subplot(4,3,[3,6]); hold on;
for i = 1:size(AREA0,1)
    XX = dataStr.(field_name).xps_hv;
    plot(XX, dataStr.(field_name).AREA0(i,:), 'kx-', 'color', pp.col.fit{i}, 'markerfacecolor', pp.col.fit{i}, 'linewidth', 1.5,'markersize', 10);
end
% -- Formatting the figure
gca_props(); % grid on;
axis([...
    min(XX)-1, max(XX)+1,...
    0, 1]);
xlabel('$$ \bf  hv [eV] $$', 'Interpreter', 'latex');
ylabel('$$ \bf  AREA0  $$', 'Interpreter', 'latex');

% Plotting BE vs measurement number
subplot(4,3,[9]); hold on;
for i = 1:size(BE,1)
    XX = dataStr.(field_name).xps_hv;
    plot(XX, dataStr.(field_name).BE(i,:), 'kx-', 'color', pp.col.fit{i}, 'markerfacecolor', pp.col.fit{i}, 'linewidth', 1.5,'markersize', 10);
    errorbar(mean(XX(:)), dataStr.(field_name).BE_mu(i), dataStr.(field_name).BE_3sig(i), dataStr.(field_name).BE_3sig(i), 0, 0, 'ks-',...
        'color', [0 0 0], 'markerfacecolor', pp.col.fit{i}, 'markersize', 10, 'linewidth', 1.5);
    BE_str = string(round(dataStr.(field_name).BE_mu(i),3)) + " (" + string(round(dataStr.(field_name).BE_3sig(i),3)) + ") eV";
    text(max(XX),dataStr.(field_name).BE_mu(i),BE_str,'Interpreter',"latex",'FontSize',15,'FontWeight','bold');
end
% -- Formatting the figure
gca_props(); % grid on;
axis([...
    min(XX)-1, max(XX)+1,...
    min(dataStr.(field_name).BE_mu(:)) - 1.25*max(dataStr.(field_name).BE_3sig(:))-0.005,...
    max(dataStr.(field_name).BE_mu(:)) + 1.25*max(dataStr.(field_name).BE_3sig(:))+0.005]);
xlabel('$$ \bf  hv [eV] $$', 'Interpreter', 'latex');
ylabel('$$ \bf  BE [eV]  $$', 'Interpreter', 'latex');

% Plotting CHISQ vs measurement number
subplot(4,3,[12]); hold on;
XX = dataStr.(field_name).xps_hv;
plot(XX, dataStr.(field_name).CHISQ, 'kx-', 'color', pp.col.fit{1}, 'markerfacecolor', pp.col.fit{1}, 'linewidth', 1.5,'markersize', 10);
errorbar(mean(XX(:)), dataStr.(field_name).CHISQ_mu, dataStr.(field_name).CHISQ_3sig, dataStr.(field_name).CHISQ_3sig, 0, 0, 'ks-',...
    'color', [0 0 0], 'markerfacecolor', pp.col.fit{1}, 'markersize', 10, 'linewidth', 1.5);
% -- Formatting the figure
gca_props(); % grid on;
axis([...
    min(XX)-1, max(XX)+1,...
    dataStr.(field_name).CHISQ_mu - 1.25*dataStr.(field_name).CHISQ_3sig-0.005,...
    dataStr.(field_name).CHISQ_mu + 1.25*dataStr.(field_name).CHISQ_3sig+0.005]);
xlabel('$$ \bf  hv [eV] $$', 'Interpreter', 'latex');
ylabel('$$ \bf  CHISQ  $$', 'Interpreter', 'latex');
CHISQ_str = string(round(dataStr.(field_name).CHISQ_mu,3)) + " (" + string(round(dataStr.(field_name).CHISQ_3sig,3)) + ")";
text(0.05,0.8,CHISQ_str,'Interpreter',"latex",'Units',"normalized",'FontSize',15,'FontWeight','bold');
% Saving the figure
print(path_save + dataStr.(field_name).filename,'-dpng', '-r500');