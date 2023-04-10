function rs_arpes_data = resample_arpes_data(arpes_data, rs_num, rs_type)
% rs_arpes_data = resample_arpes_data(arpes_data, rs_num, rs_type)
%   This is a function that resamples a collection of ARPES data using a
%   Bayesian Bootstrap method. For N independent ARPES data, the data
%   is resampled by using a linear combination of the available data each 
%   of which are scaled by coefficients (between 0 - 1) that are determined 
%   by a uniform random number generator. This sufficiently mixes all of 
%   the experimental data, to create a new 'pseudo' data set that can be
%   used for statistical analysis.
%
%   REQ. FUNCTIONS: (none)
%
%   IN:
%   -   arpes_data:    	{1×N} cell-array that contains a collection of all ARPES data, loaded and processed using ARPESGUI.
%   -   rs_num:         scalar number of the total number of resampled .
%   -   rs_type:        string of the data resampling technique to use (if empty, Default: "BayesianBootstrap" ("BB"), or "JackKnife" ("JK"))
%
%   OUT: (none, only the file is saved)
%   -   rs_arpes_data:  {1×N} cell-array that contains all the resampled data.

%% Default parameters
if nargin < 3; rs_type = "BB"; end
if nargin < 2; rs_num = 1e3; end
if isempty(rs_type); rs_type = "BB"; end
if isempty(rs_num); rs_num = 1e3; end

%% - 1 - Verifying consistency of the ARPES data grid size
% -- Extracting the ARPES data dimensions
for i = 1:length(arpes_data)
    eb_size(i)  = size(arpes_data{i}.eb, 1);
    kx_size(i)  = size(arpes_data{i}.eb, 2);
end
% -- Storing the minimum index values, to ensure all original ARPES data has a consistent domain
eb_indx_min = min(eb_size(:));  % find the minimum eb size
kx_indx_min = min(kx_size(:));  % find the minimmum kx size
% -- Cropping the original ARPES data to have a consistent domain (if required)
for i = 1:length(arpes_data)
    EB{i}   = arpes_data{i}.eb(1:eb_indx_min, 1:kx_indx_min);
    KX{i}   = arpes_data{i}.kx(1:eb_indx_min, 1:kx_indx_min);
    DATA{i}	= arpes_data{i}.data(1:eb_indx_min, 1:kx_indx_min);
end

%% - 2 - Defining the Eb and Kx dimensions for the resampled data
rs_arpes_data.Type   = "Eb(k) - Resampled";
rs_arpes_data.kx     = KX{1};
rs_arpes_data.eb     = EB{1};
    
%% - 3.0 - Creating an artificial data set based on all the ARPES data
rs_kx 	= rs_arpes_data.kx;
rs_eb   = rs_arpes_data.eb;
%% - 3.1 - Bayesian Bootstrap method
if lower(rs_type) == "bb" || lower(rs_type) == "bayesianbootstrap"
    for n = 1:rs_num
        % -- Initialising the resampled ARPES data arrays
        rs_data = [];
        rs_sf   = [];
        for i = 1:length(DATA)
            rs_sf(i)            = rand(1);
            rs_data(:,:,i)      = rs_sf(i).* DATA{i};
        end
        % --- Applying cross-correlation analysis
        [rs_data, ~]            = SumScanXC(rs_eb(:,1), rs_data, 0.5); 
        rs_data(isnan(rs_data)) = 0;
        % --- Renormalising the resampled ARPES data
        rs_data                 = rs_data - min(rs_data(:));
        rs_data                 = rs_data ./ max(rs_data(:));
        % --- Storing the ARPES data into the resampled cell-array
        rs_arpes_data.sf{n}     = rs_sf;
        rs_arpes_data.data{n}   = rs_data;
    end
elseif lower(rs_type) == "jk" || lower(rs_type) == "jackknife"
    for n = 1:length(DATA)
        % -- Initialising the resampled ARPES data arrays
        rs_data     = [];
        rs_sf       = ones(1,length(DATA));
        rs_sf(n)    = 0;
        for i = 1:length(DATA)
            rs_data(:,:,i)      = rs_sf(i).* DATA{i};
        end
        % --- Applying cross-correlation analysis
        [rs_data, ~]            = SumScanXC(rs_eb(:,1), rs_data, 0.5); 
        rs_data(isnan(rs_data)) = 0;
        % --- Renormalising the resampled ARPES data
        rs_data                 = rs_data - min(rs_data(:));
        rs_data                 = rs_data ./ max(rs_data(:));
        % --- Storing the ARPES data into the resampled cell-array
        rs_arpes_data.sf{n}     = rs_sf;
        rs_arpes_data.data{n}   = rs_data;
    end
end
end