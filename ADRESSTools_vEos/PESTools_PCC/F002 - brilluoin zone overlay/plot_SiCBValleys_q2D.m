function [] = plot_SiCBValleys_q2D(eb_cut, origin, cb_shifts)
% [] = plot_SiCBValleys_q2D(eb_cut, origin, cb_shifts)
%   This function plots the iso-energetic surfaces of elliptical valleys
%   for Silicon. The 6 elliptical valleys are now broken in degeneracy, and
%   colored differently to distinguish the 6XYZ -> 4XY + 2Z valley
%   polarisation.
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
%% 1.2 - Extracting the CB as cylinders
k0 = 0.85;       % Defining the position of the CBM in k-space
% -- IsoE for X-valleys (cylinders)
% --- Extracting unit cylindrical valleys with radius defined by ellipse
[pX_xx, pX_yy, pX_zz] = cylinder(k_calc(eb_cut,ml,y_ebshift), 1e3);
[mX_xx, mX_yy, mX_zz] = cylinder(k_calc(eb_cut,ml,y_ebshift), 1e3);
% --- Scaling the cylinder in z to be the same size as ellipses in k-space
pX_yy = pX_yy /2;
mX_yy = mX_yy /2;
pX_zz = 2*(k_calc(eb_cut,mt,z_ebshift))*(pX_zz-0.5);
mX_zz = 2*(k_calc(eb_cut,mt,z_ebshift))*(mX_zz-0.5);
% --- Shifting the valleys to the correct locations in k-space
pX_xx = pX_xx + k0 + origin(1);
pX_yy = pX_yy + 0 + origin(2);
pX_zz = pX_zz + 0 + origin(3); 
mX_xx = mX_xx - k0 + origin(1);
mX_yy = mX_yy + 0 + origin(2);
mX_zz = mX_zz + 0 + origin(3);
% -- IsoE for Y-valleys (cylinders)
% --- Extracting unit cylindrical valleys with radius defined by ellipse
[pY_xx, pY_yy, pY_zz] = cylinder(k_calc(eb_cut,ml,y_ebshift), 1e3);
[mY_xx, mY_yy, mY_zz] = cylinder(k_calc(eb_cut,ml,y_ebshift), 1e3);
% --- Scaling the cylinder in z to be the same size as ellipses in k-space
pY_xx = pY_xx /2;
mY_xx = mY_xx /2;
pY_zz = 2*(k_calc(eb_cut,mt,z_ebshift))*(pY_zz-0.5);
mY_zz = 2*(k_calc(eb_cut,mt,z_ebshift))*(mY_zz-0.5);
% --- Shifting the valleys to the correct locations in k-space
pY_xx = pY_xx + 0+origin(1);
pY_yy = pY_yy + k0+origin(2);
pY_zz = pY_zz + 0 + origin(3); 
mY_xx = mY_xx + 0+origin(1);
mY_yy = mY_yy - k0 + origin(2);
mY_zz = mY_zz + 0 + origin(3);
% -- IsoE for Z-valleys (cylinders)
% --- Extracting unit cylindrical valleys with radius defined by ellipse
[pZ_xx, pZ_yy, pZ_zz] = cylinder(k_calc(eb_cut,mt,z_ebshift), 1e3);
[mZ_xx, mZ_yy, mZ_zz] = cylinder(k_calc(eb_cut,mt,z_ebshift), 1e3);
% --- Scaling the cylinder in z to be the same size as ellipses in k-space
pZ_zz = 2*(k_calc(eb_cut,ml,z_ebshift))*(pZ_zz-0.5);
mZ_zz = 2*(k_calc(eb_cut,ml,z_ebshift))*(mZ_zz-0.5);
% --- Shifting the valleys to the correct locations in k-space
pZ_xx = pZ_xx + 0+origin(1);
pZ_yy = pZ_yy + 0+origin(2);
pZ_zz = pZ_zz + k0 + origin(3); 
mZ_xx = mZ_xx + 0+origin(1);
mZ_yy = mZ_yy + 0+origin(2);
mZ_zz = mZ_zz - k0 + origin(3);

%% - 2.0 - Plotting the isoe surface
hold on;
% -- Defining plot constants
if cb_shifts(1)==0 && cb_shifts(2)==0 && cb_shifts(3)==0
    xy_valley_col = [0.1,0.6,0.1];
    pz_valley_col = [0.6,0.1,0.6];
    mz_valley_col = pz_valley_col;
else
    xy_valley_col = [0.1,0.6,0.1];
    pz_valley_col = [0.6,0.1,0.6];
    mz_valley_col = pz_valley_col;
end
% -- Plotting X-valleys
surf(pX_xx, pX_yy, pX_zz, 'Edgecolor', 'none', 'facecolor', xy_valley_col);
surf(mX_xx, mX_yy, mX_zz, 'Edgecolor', 'none', 'facecolor', xy_valley_col);
% -- Plotting Y-valleys
surf(pY_xx, pY_yy, pY_zz, 'Edgecolor', 'none', 'facecolor', xy_valley_col);
surf(mY_xx, mY_yy, mY_zz, 'Edgecolor', 'none', 'facecolor', xy_valley_col);
% -- Plotting Z-valleys
surf(pZ_xx, pZ_yy, pZ_zz, 'Edgecolor', 'none', 'facecolor', pz_valley_col);
surf(mZ_xx, mZ_yy, mZ_zz, 'Edgecolor', 'none', 'facecolor', mz_valley_col);
% - 2.1 - Plotting the caps to the cylinders
% - X-valleys
fill3(pX_xx(1,:), pX_yy(1,:), pX_zz(1,:), xy_valley_col, 'Edgecolor', 'none');
fill3(pX_xx(end,:), pX_yy(end,:), pX_zz(end,:), xy_valley_col, 'Edgecolor', 'none');
fill3(mX_xx(1,:), mX_yy(1,:), mX_zz(1,:), xy_valley_col, 'Edgecolor', 'none');
fill3(mX_xx(end,:), mX_yy(end,:), mX_zz(end,:), xy_valley_col, 'Edgecolor', 'none');
% - Y-valleys
fill3(pY_xx(1,:), pY_yy(1,:), pY_zz(1,:), xy_valley_col, 'Edgecolor', 'none');
fill3(pY_xx(end,:), pY_yy(end,:), pY_zz(end,:), xy_valley_col, 'Edgecolor', 'none');
fill3(mY_xx(1,:), mY_yy(1,:), mY_zz(1,:), xy_valley_col, 'Edgecolor', 'none');
fill3(mY_xx(end,:), mY_yy(end,:), mY_zz(end,:), xy_valley_col, 'Edgecolor', 'none');
% - Z-valleys
fill3(pZ_xx(1,:), pZ_yy(1,:), pZ_zz(1,:), pz_valley_col, 'Edgecolor', 'none');
fill3(pZ_xx(end,:), pZ_yy(end,:), pZ_zz(end,:), pz_valley_col, 'Edgecolor', 'none');
fill3(mZ_xx(1,:), mZ_yy(1,:), mZ_zz(1,:), pz_valley_col, 'Edgecolor', 'none');
fill3(mZ_xx(end,:), mZ_yy(end,:), mZ_zz(end,:), pz_valley_col, 'Edgecolor', 'none');
% - 3.0 Formatting the figure
axis tight; axis equal; 
end