function pes_model = nlayer_pes_model(pes_hv, lyr_mat, lyr_cls, lyr_thick, theta, gamma, plot_result)
% pes_model = nlayer_pes_model(pes_hv, lyr_mat, lyr_cls, lyr_thick, theta, gamma, plot_result)
%   Function that determines the total photoelectron intensity originating
%   from N independent layers in a sample that consists of several
%   different materials of a user-defined thickness. The core-level being
%   probed is also defined, as well as the experimental geometry. 
%   --------------------------------------------------------------------
%   The model does the following:
%
%       (1)  Using the input photon energies, the electron kinetic energies
%       are determined so that the attenuation length / inelastic-mean free
%       path (IMFP) of electrons in each layer can be determined. 
%       Using:
%           Ek = Ehv - Ebe - WF ,
%       Ek is the electron kinetic energy; Ehv is the input photon energy;
%       Ebe is the binding energy of the electron, which is determined by
%       in the Photoionisation Energy and Fluorescence Database (PIEFD)
%       whose value depends on the input core-level defined; WF is the sample
%       work function, which is determined by the Materials Properties
%       Database (MPD), whose value depends on the layer material.
%
%       (2) The total intensity of photoelectrons from each layer depends
%       on several parameters, all of which can be determined;
%           I = N * F
%       where N is the number of atoms per unit volume, determined by:
%               N = N_AVOGADRO * DENSITY / MOLECULAR WEIGHT.
%       The density and molecular weight are both contained within the MPD
%       for each material defined.
%       F is the relative sensitivity factor, determined by:
%               F = XSECT * ANG_CORR * IMFP * TRANS_CORR * exp(-C/IMFP)
%       where XSECT is the photoionisation cross-section, determined within 
%       the Photoionisation Cross-Section and Asymmetry Database (PIXSAD).
%       ANG_CORR is the angular correction factor (function of the
%       ASYmmetry parameter and the ANGle between the incident x-rays and
%       analyser), which is determined by;
%                   ANG_CORR = 1 + 0.5*ASY (1.5 * sin(ANG)^2 - 1).
%       The IMFP is determined from the Tanuma et. al. TPP-2M formalism,
%       defined within the 'imfp_tpp2m_mpd()' function, only the kinetic
%       energy and material are required, which are known from the user
%       inputs.
%       TRANS_CORR is the transmission correction (function of the kinetic
%       energy) and this is assumed to be 1 for this model.
%       C is the correction for surface contamination, which is
%       proportional to layer thickness, which we assume is equal to 0 for
%       this model.
%       Thus, the TOTAL intensity I0 of each layer can be determined using
%       the formalism discussed here.
%
%       (3) The Beer-Lambert law is used to model the attenuation of the
%       photoelectron intensity as a function of depth. This model states
%       that the photoelectron intensity decreases exponentially due to 
%       inelastic scattering, with a decay constant equal to the IMFP.
%       Thus, for a bulk, single-layer system, the Beer-Lambert law states:
%           I(z) = I0 * exp(-z / IMFP),
%       which gives the photoelectron intensity, I, that emerges from a
%       depth, z, within the sample. This is applied to each layer in the
%       n-layered stack defined by the user, where I0 is determined from
%       part (2), discussed above. The total integral intensity for each
%       layer in the sample stack then gives the total emitted
%       photoelectron intensity.
%           
%       (4) Once the integral intensites from each layer are determined,
%       the final intensities are expressed in terms of a 'relative
%       contribution'. Here, each layer intensity is normalised by the sum
%       of all the intensities at each photon energy.
%   --------------------------------------------------------------------
%
%   IN:
%   -   pes_hv:     Nx1 column vector of the input photon energy.
%   -   lyr_mat:  	Mx1 cell-vector of the material for each layer in the stack; e.g. "Si", "SiO2", "Al2O3"...
%   -   lyr_cls:  	Mx1 cell-vector of the core-level probed for each layer in the stack; e.g. "Si(2p)", "Si(2s)", "Al(2p)"...
%   -   lyr_thick:	Mx1 cell-vector of the thickness of each layer in the stack (in nano-metres)
%   -   theta:      scalar of the angle between the emitted photoelectrons and the surface normal (degrees).
%   -   gamma:      scalar of the angle between the incident x-rays and the analyser (degrees).
%
%   OUT:
%   -   pes_model:      data structure that contains all the pes model parameters and variables.

%% Default parameters
% -- Defining the default parameters
if nargin < 7; plot_result = 1;  end
if nargin < 6; gamma = 70;  end
if nargin < 5; gamma = 70; theta = 0;  end
if isempty(plot_result); plot_result = 1; end
if isempty(gamma); gamma = 80; end
if isempty(theta); theta = 0; end
if isempty(pes_hv); pes_hv = linspace(350, 1500, 1e4); end
% -- Extracting the total number of layers to be probed
Nlyrs       = length(lyr_mat);
% -- Verify that the input layer variables are consistent i nsize
lyr_cls   	= lyr_cls(1:Nlyrs);
lyr_thick   = lyr_thick(1:Nlyrs);

%% 1 - Determination of the IMFP at all kinetic energy for each layer material
%% 1.1 - Defining the notation for core levels
% --- Defining in terms of orbital notation
cls_00 = {...
    '1s',...
    '2s', '2p',...
    '3s', '3p', '3d',...
    '4s', '4p', '4d', '4f',...
    '5s', '5p', '5d', '5f',...
    '6s', '6p', '6d',...
    '7s',...
    };
% --- Defining in terms of IUPAC notation
cls_01 = {...
    'K',...                             % 1s
    'L1', 'L3',...                      % 2s, 2p
    'M1', 'M3', 'M5',...                % 3s, 3p, 3d
    'N1', 'N3', 'N5', 'N7',...          % 4s, 4p, 4d, 4f
    'O1', 'O3', 'O5', 'O7',...          % 5s, 5p, 5d, 5f
    'P1', 'P3', 'P5',...                % 6s, 6p, 6d
    'Q1',...                            % 7s
    };
%% 1.2 - Extracting the IMFP
for i = 1:Nlyrs
    %% 1.2.1 - Extracting the binding energy for each layer material core-level
    % (A) Determine the core-level being probed
    lyr_elem{i}     = lyr_cls{i}(1:end-4);              % Corresponding element
    lyr_elem_cls{i} = lyr_cls{i}(end-2:end-1);          % Corresponding core-level
    lyr_piefd{i}    = get_piefd_props(lyr_elem{i});  	% Extract properties from PIEFD
    findx           = find(contains(cls_00, string(lyr_elem_cls{i})));   % Find index of the core-level
    if isempty(findx); findx = find(contains(cls_01, string(lyr_elem_cls{i}))); end
    % (B) Extracting the core-level binding energy
    cls_iupac       = cls_01(findx);                    % Core-level in IUPAC notation
    piefd_varnames  = lyr_piefd{i}.Properties.VariableNames;
    findx           = strcmpi(piefd_varnames(:), "BE_" + string(cls_iupac));
    lyr_be{i}       = lyr_piefd{i}{1,findx};            % Extracting the binding energy of the core-level
    lyr_Z{i}        = lyr_piefd{i}.ATOM_ZNUM;           % Extracting the atomic Z number of the layer
    %% 1.2.2 - Extracting the work function for each layer material
    % (A) Extract the properties from the MPD
    lyr_mpd{i}	= get_mpd_props(lyr_mat{i});
    % (B) Extract the work function, density and atomic mass
    lyr_WF{i}       = lyr_mpd{i}.ELE_WFUNC;
    lyr_density{i}	= lyr_mpd{i}.DENSITY;
    lyr_atmass{i}   = lyr_mpd{i}.ATOM_MASS;
    % -- Validity check - if the result is NaN, assume a value of zero
    if isnan(lyr_WF{i}); lyr_WF{i} = 0; end
    %% 1.2.3 - Determine the kinetic energy of photoelectrons for each layer material
    pes_ke{i}   = pes_hv - lyr_be{i} - lyr_WF{i};
    %% 1.2.4 - Determine the IMFP at all kinetic energies for each layer material
    lyr_IMFP{i} = imfp_tpp2m_mpd(pes_ke{i}, lyr_mat{i});
    % Ensuring the imfp is a row vector
    if size(lyr_IMFP{i}, 1) > 1; lyr_IMFP{i} = lyr_IMFP{i}'; end    
    % Convert the IMFP into nm
    lyr_IMFP{i} = 0.1 * lyr_IMFP{i};
end
% - Determine the average Z number for all layers
lyr_avgZ = mean(cell2mat(lyr_Z(:)));

%% 2 - Determine the maximum intensity I0 for each layer material
for i = 1:Nlyrs
    %% 2.1 - Determine the total number of atoms per unit volume
    Nav         = 6.022e23;
    lyr_N{i}    = Nav.*(lyr_density{i} ./ lyr_atmass{i});   % number of atoms per unit volume
    %% 2.2 - Extracting the photoionisation cross section
    [~, lyr_sigma{i}, lyr_beta{i}] = get_pixsad_props(lyr_elem{i}, lyr_elem_cls{i}, pes_hv);  	% Extract properties from PIXSAD
    % -- Validity check - if the result is NaN, assume a value of zero (or 2 for asymmetry)
    if isnan(lyr_sigma{i}); lyr_sigma{i} = 0; end
    if isnan(lyr_beta{i});  lyr_beta{i} = 2; end
    %% 2.3 - Determine the adjusted asymmetry factor for solid materials
    a = 0.781; b = 0.00514; c = 0.000031;
    lyr_beta_adj{i} = lyr_beta{i} .* (a - b.*lyr_avgZ + c.*lyr_avgZ.^2);    
    %% 2.4 - Determine the angular correction factor using the asymmetry parameter
    lyr_gamma{i} = 1 + 0.5.*lyr_beta_adj{i} .* (1.5.*(sin(deg2rad(gamma))).^2 - 1);
    %% 2.5 - Determine other variables for I0 calculation
    lyr_T{i}	= 1;        % Transmission correction
    lyr_c{i}	= 0;        % Surface contamination correction
    %% 2.6 - Determine the relative sensitivity factor F
    lyr_F{i}	= lyr_sigma{i}.*lyr_gamma{i}.*lyr_IMFP{i}.*lyr_T{i}.*exp(-lyr_c{i}./lyr_IMFP{i});
    %% 2.7 - Determine the value of I0
    lyr_I0{i}   = lyr_N{i} .* lyr_F{i};
end

%% 3 - Use Beer-Lambert law to determine the total photoelectron intensity for each layer
% For one layer (Bulk)
if Nlyrs == 1
    % -- Extracting the layer intensities
    lyr_ints{1}     = lyr_I0{1} .* exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta)));
    % -- Normalising each curve in terms of its relative contribution
    lyr_norm        = lyr_ints{1};
    lyr_ints0{1}    = lyr_ints{1} ./ lyr_norm;
% For two layers (S1 -> Bulk)
elseif Nlyrs == 2
    % -- Extracting the layer intensities
    lyr_ints{1}     = lyr_I0{1} .* (1 - exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta))));
    lyr_ints{2}     = lyr_I0{2} .* exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta)));
    % -- Normalising each curve in terms of its relative contribution
    lyr_norm        = lyr_ints{1} + lyr_ints{2};
    lyr_ints0{1}    = lyr_ints{1} ./ lyr_norm;
    lyr_ints0{2}    = lyr_ints{2} ./ lyr_norm;
% For three layers (S1 -> S2 -> Bulk)
elseif Nlyrs == 3
    % -- Extracting the layer intensities
    lyr_ints{1}     = lyr_I0{1} .* (1 - exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta))));
    lyr_ints{2}     = lyr_I0{2} .* (1 - exp(-lyr_thick{2} ./ (lyr_IMFP{2} .* cos(theta)))) .* exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta)));
    lyr_ints{3}     = lyr_I0{3} .* exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta))).* exp(-lyr_thick{2} ./ (lyr_IMFP{2} .* cos(theta)));
    % -- Normalising each curve in terms of its relative contribution
    lyr_norm        = lyr_ints{1} + lyr_ints{2} + lyr_ints{3};
    lyr_ints0{1}    = lyr_ints{1} ./ lyr_norm;
    lyr_ints0{2}    = lyr_ints{2} ./ lyr_norm;
    lyr_ints0{3}    = lyr_ints{3} ./ lyr_norm;
% For four layers (S1 -> S2 -> S3 -> Bulk)
elseif Nlyrs == 4
    % -- Extracting the layer intensities
    lyr_ints{1}     = lyr_I0{1} .* (1 - exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta))));
    lyr_ints{2}     = lyr_I0{2} .* (1 - exp(-lyr_thick{2} ./ (lyr_IMFP{2} .* cos(theta)))) .* exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta)));
    lyr_ints{3}     = lyr_I0{3} .* (1 - exp(-lyr_thick{3} ./ (lyr_IMFP{3} .* cos(theta)))) .* exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta))).* exp(-lyr_thick{2} ./ (lyr_IMFP{2} .* cos(theta)));
    lyr_ints{4}     = lyr_I0{4} .* exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta))).* exp(-lyr_thick{2} ./ (lyr_IMFP{2} .* cos(theta))).* exp(-lyr_thick{3} ./ (lyr_IMFP{3} .* cos(theta)));
    % -- Normalising each curve in terms of its relative contribution
    lyr_norm        = lyr_ints{1} + lyr_ints{2} + lyr_ints{3} + lyr_ints{4};
    lyr_ints0{1}    = lyr_ints{1} ./ lyr_norm;
    lyr_ints0{2}    = lyr_ints{2} ./ lyr_norm;
    lyr_ints0{3}    = lyr_ints{3} ./ lyr_norm;
    lyr_ints0{4}    = lyr_ints{4} ./ lyr_norm;
% For five layers (S1 -> S2 -> S3 -> S4 -> Bulk)
elseif Nlyrs == 5
    % -- Extracting the layer intensities
    lyr_ints{1}     = lyr_I0{1} .* (1 - exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta))));
    lyr_ints{2}     = lyr_I0{2} .* (1 - exp(-lyr_thick{2} ./ (lyr_IMFP{2} .* cos(theta)))) .* exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta)));
    lyr_ints{3}     = lyr_I0{3} .* (1 - exp(-lyr_thick{3} ./ (lyr_IMFP{3} .* cos(theta)))) .* exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta))).* exp(-lyr_thick{2} ./ (lyr_IMFP{2} .* cos(theta)));
    lyr_ints{4}     = lyr_I0{4} .* (1 - exp(-lyr_thick{4} ./ (lyr_IMFP{4} .* cos(theta)))) .* exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta))).* exp(-lyr_thick{2} ./ (lyr_IMFP{2} .* cos(theta))).* exp(-lyr_thick{3} ./ (lyr_IMFP{3} .* cos(theta)));
    lyr_ints{5}     = lyr_I0{5} .* exp(-lyr_thick{1} ./ (lyr_IMFP{1} .* cos(theta))).* exp(-lyr_thick{2} ./ (lyr_IMFP{2} .* cos(theta))).* exp(-lyr_thick{3} ./ (lyr_IMFP{3} .* cos(theta))).* exp(-lyr_thick{4} ./ (lyr_IMFP{4} .* cos(theta)));
    % -- Normalising each curve in terms of its relative contribution
    lyr_norm        = lyr_ints{1} + lyr_ints{2} + lyr_ints{3} + lyr_ints{4} + lyr_ints{5};
    lyr_ints0{1}    = lyr_ints{1} ./ lyr_norm;
    lyr_ints0{2}    = lyr_ints{2} ./ lyr_norm;
    lyr_ints0{3}    = lyr_ints{3} ./ lyr_norm;
    lyr_ints0{4}    = lyr_ints{4} ./ lyr_norm;
    lyr_ints0{5}    = lyr_ints{5} ./ lyr_norm;
% For any higher order layer, return an error
else
    msg = 'Toooooo many layers defined; only calculates a maximum of 5 layers = 4 (overlayers) + 1 (bulk).'; 
    error(msg)       
end

%% 4 - Creating a MATLAB data structure to store all of the results
% - Defining the structure
pes_model               = struct();
% - Defining the input arguments
pes_model.pes_hv        = pes_hv;
pes_model.lyr_mat       = lyr_mat;
pes_model.lyr_cls       = lyr_cls;
pes_model.lyr_thick     = lyr_thick;
pes_model.theta         = theta;
pes_model.gamma         = gamma;
% - Defining the material properties
pes_model.Nlyrs         = Nlyrs;
pes_model.lyr_mpd       = lyr_mpd;
pes_model.lyr_piefd     = lyr_piefd;
% - Defining the properties for extracting the IMFP
pes_model.lyr_elem      = lyr_elem;
pes_model.lyr_elem_cls	= lyr_elem_cls;
pes_model.lyr_be        = lyr_be;
pes_model.lyr_Z         = lyr_Z;
pes_model.lyr_avgZ      = lyr_avgZ;
pes_model.lyr_WF        = lyr_WF;
pes_model.lyr_density 	= lyr_density;
pes_model.lyr_atmass   	= lyr_atmass;
pes_model.pes_ke        = pes_ke;
pes_model.lyr_IMFP      = lyr_IMFP;
% - Defining the properties for extracting value of I0
pes_model.lyr_N         = lyr_N;
pes_model.lyr_sigma     = lyr_sigma;
pes_model.lyr_beta      = lyr_beta;
pes_model.lyr_beta_adj	= lyr_beta_adj;
pes_model.lyr_gamma     = lyr_gamma;
pes_model.lyr_T         = lyr_T;
pes_model.lyr_c         = lyr_c;
pes_model.lyr_F         = lyr_F;
pes_model.lyr_I0        = lyr_I0;
% - Defining the final intensities
pes_model.lyr_ints    	= lyr_ints;
pes_model.lyr_ints0  	= lyr_ints0;

%% 5 - Plotting the multilayer model stack and solutions
if plot_result == 1; view_nlayer_pes_model(pes_model); end

end