function [] = plot_GeCBValleys_3D(eb_cut, origin, cb_shifts)
% [] = plot_SiCBValleys_3D(eb_cut, origin, cb_shifts)
%   This function plots the iso-energetic surfaces of elliptical valleys
%   for Silicon. This is 6 elliptical valleys, 85% along the Gamma-X high
%   symmetry line.
%
%   IN:
%   -  eb_cut:     	scalar of the energy above the CBM to slice for the iso-surface [eV].
%   -  origin:      1x3 row vector of the origin location of the Gama point [x,y,z].
%   -  cb_shifts: 	1x3 row vector of the relative energy shift of each valley [eb_x, eb_y, eb_z].
%
%   OUT: (figure)

%% Defining the default parameters
if nargin < 1; eb_cut = 0.2; origin = [0, 0, 0]; cb_shifts = [0, 0, 0]; end
if nargin < 2; origin = [0, 0, 0]; cb_shifts = [0, 0, 0];  end
if nargin < 3; cb_shifts = [0, 0, 0]; end
if isempty(eb_cut); eb_cut = 0.2;  end
if isempty(origin); origin = [0, 0, 0];  end
if isempty(cb_shifts); cb_shifts = [0, 0, 0];  end

%% 1.0 - Defining the constant parameters
x_ebshift = cb_shifts(1); y_ebshift = cb_shifts(2); z_ebshift = cb_shifts(3); 
m0 = 9.11e-31;
hbar = 4.135e-15;
ml = 1.64*m0;
mt = 0.082*m0;
%% 1.1 - Defining the dispersions of the CB parabola's
k_calc =  @(eb, eff_mass, eb_shift) sqrt(((2*eff_mass)/(hbar^2)) * (eb - eb_shift));
%% 1.2 - Extracting the CB ellipsoids
k0 = sqrt(3)*pi/(5.657);       % Defining the position of the CBM in k-space
% -- IsoE for x valleys
[X1_xx, X1_yy, X1_zz] = ellipsoid(k0 +origin(1), 0+origin(2), 0+origin(3), k_calc(eb_cut,ml,x_ebshift),k_calc(eb_cut,mt,x_ebshift),k_calc(eb_cut,mt,x_ebshift), 5e2);
[X2_xx, X2_yy, X2_zz] = ellipsoid(-k0+origin(1), 0+origin(2), 0+origin(3), k_calc(eb_cut,ml,x_ebshift),k_calc(eb_cut,mt,x_ebshift),k_calc(eb_cut,mt,x_ebshift), 5e2);
[X3_xx, X3_yy, X3_zz] = ellipsoid(origin(1), k0+origin(2), 0+origin(3), k_calc(eb_cut,mt,x_ebshift),k_calc(eb_cut,ml,x_ebshift),k_calc(eb_cut,mt,x_ebshift), 5e2);
[X4_xx, X4_yy, X4_zz] = ellipsoid(origin(1), -k0+origin(2), 0+origin(3), k_calc(eb_cut,mt,x_ebshift),k_calc(eb_cut,ml,x_ebshift),k_calc(eb_cut,mt,x_ebshift), 5e2);

%% - 2.0 - Plotting the isoe surface
hold on;
% -- Defining plot constants
valley_col = [0.8, 0.4, 0];
% -- Plotting M-valleys
% - #1
h1 = surf(X1_xx, X1_yy, X1_zz, 'Edgecolor', 'none', 'facecolor', valley_col); 
T1 = hgtransform('Parent',gca); h1.Parent=T1; T1.Matrix=makehgtform('xrotate',deg2rad(-15), 'yrotate',deg2rad(-35.264), 'zrotate',-pi/4);
% - #2
h2 = surf(X2_xx, X2_yy, X2_zz, 'Edgecolor', 'none', 'facecolor', valley_col); 
T2 = hgtransform('Parent',gca); h2.Parent=T2; T2.Matrix=makehgtform('xrotate',deg2rad(-15), 'yrotate',deg2rad(-35.264), 'zrotate',-pi/4);
% - #3
h3 = surf(X1_xx, X1_yy, X1_zz, 'Edgecolor', 'none', 'facecolor', valley_col); 
T3 = hgtransform('Parent',gca); h3.Parent=T3; T3.Matrix=makehgtform('xrotate',deg2rad(-15)+pi, 'yrotate',deg2rad(-35.264), 'zrotate',-pi/4);
% - #4
h4 = surf(X2_xx, X2_yy, X2_zz, 'Edgecolor', 'none', 'facecolor', valley_col); 
T4 = hgtransform('Parent',gca); h4.Parent=T4; T4.Matrix=makehgtform('xrotate',deg2rad(-15)+pi, 'yrotate',deg2rad(-35.264), 'zrotate',-pi/4);
% - #5
h5 = surf(X3_xx, X3_yy, X3_zz, 'Edgecolor', 'none', 'facecolor', valley_col); 
T5 = hgtransform('Parent',gca); h5.Parent=T5; T5.Matrix=makehgtform('xrotate',deg2rad(15),'yrotate',deg2rad(-35.264), 'zrotate',deg2rad(-45));
% - #6
h6 = surf(X4_xx, X4_yy, X4_zz, 'Edgecolor', 'none', 'facecolor', valley_col); 
T6 = hgtransform('Parent',gca); h6.Parent=T6; T6.Matrix=makehgtform('xrotate',deg2rad(15),'yrotate',deg2rad(-35.264), 'zrotate',deg2rad(-45));
% - #7
h7 = surf(X3_xx, X3_yy, X3_zz, 'Edgecolor', 'none', 'facecolor', valley_col); 
T7 = hgtransform('Parent',gca); h7.Parent=T7; T7.Matrix=makehgtform('xrotate',deg2rad(15)+pi,'yrotate',deg2rad(-35.264), 'zrotate',deg2rad(-45));
% - #8
h8 = surf(X4_xx, X4_yy, X4_zz, 'Edgecolor', 'none', 'facecolor', valley_col); 
T8 = hgtransform('Parent',gca); h8.Parent=T8; T8.Matrix=makehgtform('xrotate',deg2rad(15)+pi,'yrotate',deg2rad(-35.264), 'zrotate',deg2rad(-45));
% - 3.0 Formatting the figure
axis tight; axis equal;
end