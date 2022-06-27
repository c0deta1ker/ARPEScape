function [] = plot_fccBZ_3D(a_fcc, origin)
% [] = plot_fccBZ_3D(a_fcc, origin)
%   This is a function that plots the Brilluoin zone of an FCC
%   lattice, with a lattice constant a, origin location at [x0, y0, z0].
%
%   IN:
%   -  a_fcc:     	scalar of the lattice constant of the FCC lattice [nm].
%   -  origin:      1x3 row vector of the origin location of the BZ [x,y,z]
%
%   OUT: (figure)

%% Setting the default parameters
if nargin < 1; a_fcc = 5.43; origin = [0 0 0]; end
if nargin < 2; origin = [0 0 0]; end
if isempty(a_fcc); a_fcc = 5.43;  end
if isempty(origin); origin = [0 0 0];  end

%% Defining the functions constants
% Defining the brilluoin zone plot constants
lineCol = [0 0 0];
lineStyle = '-';
lineWidth = 2;
% Defining the variable constants
a_reci = 2*pi / a_fcc;

%% 1.0 - Defining the polygons of the truncated octahedron
%-Defining the hexagons in terms of a
hex = [
    [1, 0.5, 0, 0, 0.5, 1, 1];
    [0.5, 1, 1, 0.5, 0, 0, 0.5];
    [0, 0, 0.5, 1, 1, 0.5, 0]
    ]*a_reci;
%-Defining the squares in terms of a
squ = [
    [1, 1, 1, 1, 1];
    [0.5, 0, -0.5, 0, 0.5];
    [0, 0.5, 0, -0.5, 0]
    ]*a_reci;

%% 1.1 - Symmetry relationships of all the polygons
%-Extracting all of the permutations of the hexagons
h1 = [hex(1,:)+origin(1); hex(2,:)+origin(2); hex(3,:)+origin(3)];
h2 = [hex(1,:)+origin(1); -1*hex(2,:)+origin(2); hex(3,:)+origin(3)];
h3 = [-1*hex(1,:)+origin(1); hex(2,:)+origin(2); hex(3,:)+origin(3)];
h4 = [-1*hex(1,:)+origin(1); -1*hex(2,:)+origin(2); hex(3,:)+origin(3)];
h5 = [hex(1,:)+origin(1); hex(2,:)+origin(2); -1*hex(3,:)+origin(3)];
h6 = [hex(1,:)+origin(1); -1*hex(2,:)+origin(2); -1*hex(3,:)+origin(3)];
h7 = [-1*hex(1,:)+origin(1); hex(2,:)+origin(2); -1*hex(3,:)+origin(3)];
h8 = [-1*hex(1,:)+origin(1); -1*hex(2,:)+origin(2); -1*hex(3,:)+origin(3)];
%-Extracting all of the permutations of the squares
s1 = [squ(1,:)+origin(1); squ(2,:)+origin(2); squ(3,:)+origin(3)];
s2 = [-1*squ(1,:)+origin(1); squ(2,:)+origin(2); squ(3,:)+origin(3)];
s3 = [squ(2,:)+origin(1); squ(1,:)+origin(2); -1*squ(3,:)+origin(3)];
s4 = [squ(2,:)+origin(1); -1*squ(1,:)+origin(2); squ(3,:)+origin(3)];
s5 = [squ(2,:)+origin(1); squ(3,:)+origin(2); squ(1,:)+origin(3)];
s6 = [squ(2,:)+origin(1); squ(3,:)+origin(2); -1*squ(1,:)+origin(3)];
%-Extracting the faces as a convex hull
x = [h1(1,:), h2(1,:), h3(1,:), h4(1,:), h5(1,:), h6(1,:), h7(1,:), h8(1,:)];
y = [h1(2,:), h2(2,:), h3(2,:), h4(2,:), h5(2,:), h6(2,:), h7(2,:), h8(2,:)];
z = [h1(3,:), h2(3,:), h3(3,:), h4(3,:), h5(3,:), h6(3,:), h7(3,:), h8(3,:)];

%% 2.0 - Plotting the truncated octahedron Brilluoin zone
hold on;
%--Plotting the hexagons
plot3(h1(1,:), h1(2,:), h1(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
plot3(h2(1,:), h2(2,:), h2(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
plot3(h3(1,:), h3(2,:), h3(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
plot3(h4(1,:), h4(2,:), h4(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
plot3(h5(1,:), h5(2,:), h5(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
plot3(h6(1,:), h6(2,:), h6(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
plot3(h7(1,:), h7(2,:), h7(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
plot3(h8(1,:), h8(2,:), h8(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
%--Plotting the squares
plot3(s1(1,:), s1(2,:), s1(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
plot3(s2(1,:), s2(2,:), s2(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
plot3(s3(1,:), s3(2,:), s3(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
plot3(s4(1,:), s4(2,:), s4(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
plot3(s5(1,:), s5(2,:), s5(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
plot3(s6(1,:), s6(2,:), s6(3,:), 'linewidth', lineWidth,...
    'linestyle', lineStyle, 'color', lineCol);
end