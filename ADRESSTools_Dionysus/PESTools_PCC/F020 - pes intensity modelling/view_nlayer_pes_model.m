function fig = view_nlayer_pes_model(pes_model, lyr_cols)
% fig = view_nlayer_pes_model(pes_model)
%   This function plots the solutions to the 'nlayer_pes_model()' function;
%   the n-layered sample stack and the corresponding photoelectron
%   contribution from each one of the layers. The user can input the colors
%   for each one of the layers using the 'lyr_cols' argument; this is
%   convenient when you want to color-match the schematic to your own
%   drawings.
%
%   IN:
%   -   pes_model:      data structure that contains all the pes model parameters and variables (from 'nlayer_pes_model()').
%   -   lyr_cols:       Mx1 cell-vector of the [R,G,B] color of each independent layer.
%
%   OUT:
%   -   fig:            MATLAB figure object with the ARPES data plotted.

%% Default parameters
% -- Extracting plot properties
pp     	= plot_props();
% -- Defining the default parameters
if nargin < 2; lyr_cols = pp.col.fit; end
if isempty(lyr_cols); lyr_cols = pp.col.fit; end

%% 1 - Plotting the the model solutions
fig = figure();
fig.Position(3) = 2.50*500; 
fig.Position(4) = 1.0*400;
%% 1.1 - Plotting the intensity profiles for each layer
subplot(121); hold on;
gca_props(0);
x_width	= 5;
if pes_model.Nlyrs == 1
    y_cum   = cell2mat(pes_model.lyr_thick);
else
    y_cum   = cumsum(cell2mat(pes_model.lyr_thick));
    y_cum(isinf(y_cum)) = y_cum(end-1) + 10;
end
% -- Plotting each layer from bottom-up
for i = pes_model.Nlyrs:-1:1
    patch([-1, -1, 1, 1, -1].*x_width, [0, y_cum(i), y_cum(i), 0, 0].*-1,...
        lyr_cols{i}, 'edgecolor', [0 0 0]);
end
% -- Box Styling and Axis properties
ylabel('$$ \bf  Depth\ From\ Surface\ [nm] $$', 'Interpreter', 'latex');
title('Sample model stack');
axis square; 
axis([-x_width, x_width, -1*max(y_cum(:)), 0]);
% -- Adding text for each material type
for i = 1:pes_model.Nlyrs
    if i == 1;  y_loc = 0 - 0.5*(y_cum(i) - 0);
    else;       y_loc = -y_cum(i-1) - 0.5*(y_cum(i) - y_cum(i-1));
    end
    if i == pes_model.Nlyrs
        text(1.02*x_width, y_loc, "(B)"+string(pes_model.lyr_mat{i})+"["+string(pes_model.lyr_cls{i})+"]",...
        'color', 'k', 'horizontalalignment', 'left', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',11);
    else
        text(1.02*x_width, y_loc, sprintf("(S%i)",i)+string(pes_model.lyr_mat{i})+"["+string(pes_model.lyr_cls{i})+"]",...
        'color', 'k', 'horizontalalignment', 'left', 'verticalalignment', 'middle', 'FontWeight','bold', 'FontSize',11);
    end
end
%% 4.2 - Plotting the intensity profiles for each layer
subplot(122); hold on;
for i = 1:pes_model.Nlyrs
    plot(pes_model.pes_hv, pes_model.lyr_ints0{i}, 'k-', 'color', lyr_cols{i}, 'markerfacecolor', lyr_cols{i}, 'linewidth', 3.50);
    max_vals(i) = max(pes_model.lyr_ints0{i}(:));
end
legend(pes_model.lyr_mat, 'location', 'best');
gca_props(0);
title('Photoelectron intensity contributions');
xlabel('$$ \bf  Photon\ energy\ [eV] $$', 'Interpreter', 'latex');
ylabel('$$ \bf  Relative\ Contribution $$', 'Interpreter', 'latex');
ylim([0, 1.05.*max(max_vals(:))]);
ax = gca;
ax.YAxisLocation = 'right';


end