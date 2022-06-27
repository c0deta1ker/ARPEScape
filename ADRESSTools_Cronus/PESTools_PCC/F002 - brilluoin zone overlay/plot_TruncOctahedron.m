function [] = plot_TruncOctahedron(a_fcc, origin, colTrunc, colPyra)
% [] = plot_TruncOctahedron(a_fcc, origin, colTrunc, colPyra)
%   This is a function that plots a construction of the
%   truncated octahedron, showing the regular octahedron
%   and the six right-square pyramids that are cut from it.
%
%   IN:
%   -  a_fcc:     	scalar of the lattice constant of the FCC lattice [nm].
%   -  origin:      1x3 row vector of the origin location of the truncated octahedron shape [x,y,z].
%   -  colTrunc: 	1x3 row vector of the color for the truncated part.
%   -  colPyra: 	1x3 row vector of the color for the pyramid caps part.
%
%   OUT: (figure)

%% Setting the default input parameters
if nargin < 1; a_fcc = 5.43; origin = [0 0 0]; colTrunc = [0.2 0.2 0.6]; colPyra = [0.6 0.2 0.2]; end
if nargin < 2; origin = [0 0 0]; colTrunc = [0.2 0.2 0.6]; colPyra = [0.6 0.2 0.2]; end
if nargin < 3; colTrunc = [0.2 0.2 0.6]; colPyra = [0.6 0.2 0.2]; end
if nargin < 4; colPyra = [0.6 0.2 0.2]; end
if isempty(a_fcc); a_fcc = 5.43;  end
if isempty(origin); origin = [0 0 0];  end
if isempty(colTrunc); colTrunc = [0.2 0.2 0.6];  end
if isempty(colPyra); colPyra = [0.6 0.2 0.2];  end

%% Defining the functions constants
%Defining the variable constants
a_reci = 2*pi / a_fcc;
%Defining the plotting constants
fontSize = 20;
lineWidth = 2.;

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

%% 1.1 - Symmetry relationships of all the truncated octahedron
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
x_oct = [h1(1,:), h2(1,:), h3(1,:), h4(1,:), h5(1,:), h6(1,:), h7(1,:), h8(1,:)];
y_oct = [h1(2,:), h2(2,:), h3(2,:), h4(2,:), h5(2,:), h6(2,:), h7(2,:), h8(2,:)];
z_oct = [h1(3,:), h2(3,:), h3(3,:), h4(3,:), h5(3,:), h6(3,:), h7(3,:), h8(3,:)];
faces_oct = convhull(x_oct, y_oct, z_oct);

%% 2.0 - Defining the polygons of the right-sqaure pyramids
%-Defining the triangles in terms of a
tri = [
    [1,     1,       1,         1,      1,      1.5,      1,        1,      1.5,    1];
    [0,     0.5,    0,      -0.5,   0,      0,          0,      0.5,    0,      -0.5];
    [0.5,  0,       -0.5,   0,      0.5,   0,          -0.5,   0,      0,      0]
    ]* a_reci;

%% 2.1 - Symmetry relationships of all the right-square pyramids
%-Extracting all of the permutations of the pyramids
p1 = [tri(1,:)+origin(1); tri(2,:)+origin(2); tri(3,:)+origin(3)];
p2 = [-1*tri(1,:)+origin(1); tri(2,:)+origin(2); tri(3,:)+origin(3)];
p3 = [tri(2,:)+origin(1); tri(1,:)+origin(2); -1*tri(3,:)+origin(3)];
p4 = [tri(2,:)+origin(1); -1*tri(1,:)+origin(2); tri(3,:)+origin(3)];
p5 = [tri(2,:)+origin(1); tri(3,:)+origin(2); tri(1,:)+origin(3)];
p6 = [tri(2,:)+origin(1); tri(3,:)+origin(2); -1*tri(1,:)+origin(3)];
%-Extracting the faces as a convex hull
x_pyr = [p1(1,:), p2(1,:), p3(1,:), p4(1,:), p5(1,:), p6(1,:)];
y_pyr = [p1(2,:), p2(2,:), p3(2,:), p4(2,:), p5(2,:), p6(2,:)];
z_pyr = [p1(3,:), p2(3,:), p3(3,:), p4(3,:), p5(3,:), p6(3,:)];
faces_pyr = convhull(x_pyr, y_pyr, z_pyr);


%% 3.0 - Plotting the truncated octahedron Brilluoin zone
hold on;
%-Plotting the 3D Brilluoin zone shades
trisurf(faces_pyr,x_pyr,y_pyr,z_pyr,'FaceColor',colPyra,'EdgeColor','none')
trisurf(faces_oct,x_oct,y_oct,z_oct,'FaceColor',colTrunc,'EdgeColor','none')

%--Plotting the hexagons
plot3(h1(1,:), h1(2,:), h1(3,:), 'k-', 'linewidth', lineWidth);
plot3(h2(1,:), h2(2,:), h2(3,:), 'k-', 'linewidth', lineWidth);
plot3(h3(1,:), h3(2,:), h3(3,:), 'k-', 'linewidth', lineWidth);
plot3(h4(1,:), h4(2,:), h4(3,:), 'k-', 'linewidth', lineWidth);
plot3(h5(1,:), h5(2,:), h5(3,:), 'k-', 'linewidth', lineWidth);
plot3(h6(1,:), h6(2,:), h6(3,:), 'k-', 'linewidth', lineWidth);
plot3(h7(1,:), h7(2,:), h7(3,:), 'k-', 'linewidth', lineWidth);
plot3(h8(1,:), h8(2,:), h8(3,:), 'k-', 'linewidth', lineWidth);
%--Plotting the squares
plot3(s1(1,:), s1(2,:), s1(3,:), 'k-', 'linewidth', lineWidth);
plot3(s2(1,:), s2(2,:), s2(3,:), 'k-', 'linewidth', lineWidth);
plot3(s3(1,:), s3(2,:), s3(3,:), 'k-', 'linewidth', lineWidth);
plot3(s4(1,:), s4(2,:), s4(3,:), 'k-', 'linewidth', lineWidth);
plot3(s5(1,:), s5(2,:), s5(3,:), 'k-', 'linewidth', lineWidth);
plot3(s6(1,:), s6(2,:), s6(3,:), 'k-', 'linewidth', lineWidth);
%--Plotting the pyramids
plot3(p1(1,:), p1(2,:), p1(3,:), 'w-', 'linewidth', lineWidth);
plot3(p2(1,:), p2(2,:), p2(3,:), 'w-', 'linewidth', lineWidth);
plot3(p3(1,:), p3(2,:), p3(3,:), 'w-', 'linewidth', lineWidth);
plot3(p4(1,:), p4(2,:), p4(3,:), 'w-', 'linewidth', lineWidth);
plot3(p5(1,:), p5(2,:), p5(3,:), 'w-', 'linewidth', lineWidth);
plot3(p6(1,:), p6(2,:), p6(3,:), 'w-', 'linewidth', lineWidth);

%-Formatting the figure
%--Adding axes labels and limits
grid on;
xlabel(' $ k_x (\AA^{-1}) $ ', 'Interpreter', 'latex');
ylabel(' $ k_y (\AA^{-1}) $ ', 'Interpreter', 'latex');
zlabel(' $ k_z (\AA^{-1}) $ ', 'Interpreter', 'latex');
axis([[-1, 1], [-1, 1], [-1, 1]]*1.0*a_reci);
set(gca,'fontsize', fontSize, 'xtick', -20:1:20, 'ytick', -20:1:20, 'ztick', -20:1:20);
%--Setting the axis and camera lighting
axis tight;  axis vis3d; axis off; camlight(-112,24);
pbaspect([1 1 1]);
view(-158, 16);
end