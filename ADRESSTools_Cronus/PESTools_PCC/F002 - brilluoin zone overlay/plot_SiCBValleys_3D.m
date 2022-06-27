function [] = plot_SiCBValleys_3D(eb_cut, origin, cb_shifts)
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
ml = 0.92*m0;
mt = 0.19*m0;
%% 1.1 - Defining the dispersions of the CB parabola's
k_calc =  @(eb, eff_mass, eb_shift) sqrt(((2*eff_mass)/(hbar^2)) * (eb - eb_shift));
%% 1.2 - Extracting the CB ellipsoids
k0 = 0.85;       % Defining the position of the CBM in k-space
% -- IsoE for x valleys
[pX_xx, pX_yy, pX_zz] = ellipsoid(k0+origin(1), 0+origin(2), 0+origin(3), k_calc(eb_cut,ml,x_ebshift),k_calc(eb_cut,mt,x_ebshift),k_calc(eb_cut,mt,x_ebshift), 1e3);
[mX_xx, mX_yy, mX_zz] = ellipsoid(-k0+origin(1), 0+origin(2), 0+origin(3), k_calc(eb_cut,ml,x_ebshift),k_calc(eb_cut,mt,x_ebshift),k_calc(eb_cut,mt,x_ebshift), 1e3);
% -- IsoE for y valleys
[pY_xx, pY_yy, pY_zz] = ellipsoid(0+origin(1),k0+origin(2),0+origin(3),k_calc(eb_cut,mt,y_ebshift),k_calc(eb_cut,ml,y_ebshift),k_calc(eb_cut,mt,y_ebshift), 1e3);
[mY_xx, mY_yy, mY_zz] = ellipsoid(0+origin(1),-k0+origin(2),0+origin(3),k_calc(eb_cut,mt,y_ebshift),k_calc(eb_cut,ml,y_ebshift),k_calc(eb_cut,mt,y_ebshift), 1e3);
% -- IsoE for z valleys
[pZ_xx, pZ_yy, pZ_zz] = ellipsoid(0+origin(1),0+origin(2),k0+origin(3),k_calc(eb_cut,mt,z_ebshift),k_calc(eb_cut,mt,z_ebshift),k_calc(eb_cut,ml,z_ebshift), 1e3);
[mZ_xx, mZ_yy, mZ_zz] = ellipsoid(0+origin(1),0+origin(2),-k0+origin(3),k_calc(eb_cut,mt,z_ebshift),k_calc(eb_cut,mt,z_ebshift),k_calc(eb_cut,ml,z_ebshift), 1e3);

%% - 2.0 - Plotting the isoe surface
hold on;
% -- Defining plot constants
if cb_shifts(1)==0 && cb_shifts(2)==0 && cb_shifts(3)==0
    xy_valley_col = [0.1,0.6,0.1];
    z_valley_col = [0.1,0.6,0.1];
else
    xy_valley_col = [0.1,0.6,0.1];
    z_valley_col = [0.6,0.1,0.6];
end
% -- Plotting X-valleys
surf(pX_xx, pX_yy, pX_zz, 'Edgecolor', 'none', 'facecolor', xy_valley_col);
surf(mX_xx, mX_yy, mX_zz, 'Edgecolor', 'none', 'facecolor', xy_valley_col);
% -- Plotting Y-valleys
surf(pY_xx, pY_yy, pY_zz, 'Edgecolor', 'none', 'facecolor', xy_valley_col);
surf(mY_xx, mY_yy, mY_zz, 'Edgecolor', 'none', 'facecolor', xy_valley_col);
% -- Plotting Z-valleys
surf(pZ_xx, pZ_yy, pZ_zz, 'Edgecolor', 'none', 'facecolor', z_valley_col);
surf(mZ_xx, mZ_yy, mZ_zz, 'Edgecolor', 'none', 'facecolor', z_valley_col);
% - 3.0 Formatting the figure
axis tight; axis equal;
end