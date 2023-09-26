close all; clear all;
path_data   = 'D:\OneDrive\PCC_MS_Project\00_SLS_Software\ARPEScape\ARPEScape-v8.1.0\PESTools_PCC\Examples\Example_Data\ADRESS_data\';
%% *(1) ARPES - 2D Data*
close all; 
% Load in ARPES data
file_names      = {...
    '001_ARPES_HSCUT.h5',...
    };
arpes_dat = cell(1,length(file_names));
for i = 1:length(file_names); arpes_dat{i} = load_adress_data(file_names{i}, path_data); end
for i = 1:length(arpes_dat); arpes_dat{i} = xcorr_scans(arpes_dat{i}); end
view_arpes_data(arpes_dat{1});
% Load in reference data
file_names_ref      = {...
    '002_ARPES_HSCUT_EF.h5',...
    };
arpes_ref = cell(1,length(file_names_ref));
for i = 1:length(file_names_ref); arpes_ref{i} = load_adress_data(file_names_ref{i}, path_data); end
for i = 1:length(arpes_ref); arpes_ref{i} = xcorr_scans(arpes_ref{i}); end
view_arpes_data(arpes_ref{1});
%% *(2) ARPES - 3D Data*
close all; 
arpes_raw = load_adress_data('003_ARPES_kxkz.h5', path_data);
arpes_dat = delete_scans(arpes_raw, {[2:2:2e3]});   % deleting every even numbered scan, so only the odd-numbered ARPES scans remain
arpes_ref = delete_scans(arpes_raw, {[1:2:2e3]});   % deleting every odd numbered scan, so only the even-numbered XPS scans remain
% -- Viewing the ARPES data
slice_args = {820, -1.0, 0.0};
view_arpes_data(arpes_dat,slice_args);
slice_args = {820, -2.7, 0.0};
view_arpes_data(arpes_ref,slice_args);