function fig = view_nlayer_sample(lyr_mat, lyr_thick, lyr_cols)
% fig = view_nlayer_sample(lyr_mat, lyr_thick, lyr_cols)
%   This function plots the solutions to the 'nlayer_pes_model()' function;
%   the n-layered sample stack and the corresponding photoelectron
%   contribution from each one of the layers. The user can input the colors
%   for each one of the layers using the 'lyr_cols' argument; this is
%   convenient when you want to color-match the schematic to your own
%   drawings.
%
%   IN:
%   -   lyr_mat:        Mx1 cell-vector of the material for each layer in the stack; e.g. "Si", "SiO2", "Al2O3"...
%   -   lyr_thick:      Mx1 cell-vector of the thickness of each layer in the stack (in nano-metres)
%   -   lyr_cols:       Mx1 cell-vector of the [R,G,B] color of each independent layer.
%
%   OUT:
%   -   fig:            MATLAB figure object with the ARPES data plotted.

%% Default parameters
% -- Extracting plot properties
pp     	= plot_props();
% -- Defining the default parameters
if nargin < 3; lyr_cols = pp.col.fit; end
if isempty(lyr_cols); lyr_cols = pp.col.fit; end
% -- Extracting the total number of layers to be probed
Nlyrs       = length(lyr_mat);

%% 1 - Plotting the the model solutions
fig = figure();
fig.Position(3) = 1.00*pp.fig5x4(1); 
fig.Position(4) = 1.0*pp.fig5x4(2);
hold on;
%% 2 - Plotting the intensity profiles for each layer
gca_props(0);
x_width	= 5;
if Nlyrs == 1
    y_cum   = cell2mat(lyr_thick);
else
    y_cum   = cumsum(cell2mat(lyr_thick));
    y_cum(isinf(y_cum)) = y_cum(end-1) + 10;
end
% -- Plotting each layer from bottom-up
for i = Nlyrs:-1:1
    patch([-1, -1, 1, 1, -1].*x_width, [0, y_cum(i), y_cum(i), 0, 0].*-1,...
        lyr_cols{i}, 'edgecolor', [0 0 0]);
end
% -- Box Styling and Axis properties
ylabel('Depth from surface (nm)');
title('Sample model stack');
axis square; 
axis([-x_width, x_width, -1*max(y_cum(:)), 0]);
% -- Adding text for each material type
for i = 1:Nlyrs
    if i == 1;  y_loc = 0 - 0.5*(y_cum(i) - 0);
    else;       y_loc = -y_cum(i-1) - 0.5*(y_cum(i) - y_cum(i-1));
    end
    if i == Nlyrs
        text(1.02*x_width, y_loc, "(B)"+string(lyr_mat{i}),...
        'fontsize', 12, 'color', 'k', 'horizontalalignment', 'left', 'verticalalignment', 'middle');  
    else
        text(1.02*x_width, y_loc, sprintf("(S%i)",i)+string(lyr_mat{i}),...
        'fontsize', 12, 'color', 'k', 'horizontalalignment', 'left', 'verticalalignment', 'middle');  
    end
end

end