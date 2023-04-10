function view_eimfp_props(material)
% view_eimfp_props(material)
%   This is a function that plots the eIMFP curves of a particular material
%   defined by the user using all availabe calculators. This can be used to
%   quickly view all the available eIMFP data for a particular element /
%   material.
%
%   IN:
%   -   material:  	string of the material whose imfp is to be shown; e.g. "Si", "SiO2", "Al2O3"...
%
%   OUT: (figure output)

%% Default parameters (Parameters for Silicon)
if nargin < 1; material = "Si";  end
if isempty(material); material = "Si"; end
%% - 1 - Extracting the material parameters from the materials database
% - Extracting the material properties
material_props = get_mpd_props(material);
% - Extracting the material properties required for the S1 formalism
rho     = material_props.DENSITY;
Nv      = material_props.ELECT_VALENCY;
M       = material_props.ATOM_MASS;
Egap    = material_props.ELE_BGAP;
Z       = material_props.ATOM_ZNUM;
%% - 2 - Extracting the eIMFP properties from each calculator
ele_imfp                 = struct();
% -- Defining the kinetic energy range
ele_imfp.Ek              = linspace(1, 5e5, 1e5);
% -- Optical data
[ele_imfp.imfp_opt, ele_imfp.dimfp_opt] = eimfp_optical(ele_imfp.Ek, material);
% -- S1 & S2 predictive equations
ele_imfp.imfp_S1         = eimfp_S1_mpd(ele_imfp.Ek, material);
ele_imfp.imfp_S2_avg     = eimfp_S2_avg(ele_imfp.Ek, Z);
% -- TPP2M predictive equations
ele_imfp.imfp_tpp2m      = eimfp_tpp2m_mpd(ele_imfp.Ek, material);
ele_imfp.imfp_tpp2m_avg  = eimfp_tpp2m_avg(ele_imfp.Ek);
% -- Unversal curve
ele_imfp.imfp_universal  = eimfp_universal(ele_imfp.Ek);
%% - 3 - Plotting the eIMFP figure
% -- Initialising the figure
fig = figure(); hold on;
fig.Position(3) = 600; 
fig.Position(4) = 450;
grid on; grid minor;
% - PLOTTING OPTICAL/EXPERIMENTAL DATA
errorbar(ele_imfp.Ek, ele_imfp.imfp_opt, ele_imfp.dimfp_opt, ele_imfp.dimfp_opt, 'rx-');
% - PLOTTING TPP2M-CURVE
loglog(ele_imfp.Ek, ele_imfp.imfp_tpp2m, ':', 'color', 'b', 'LineWidth', 2.0);
% - PLOTTING TPP2M-AVERAGE-CURVE
loglog(ele_imfp.Ek, ele_imfp.imfp_tpp2m_avg, '--', 'color', 'b', 'LineWidth', 2.0);
% - PLOTTING S1-CURVE
loglog(ele_imfp.Ek, ele_imfp.imfp_S1, ':', 'color', 'g', 'LineWidth', 2.0);
% - PLOTTING S2-AVERAGE-CURVE
loglog(ele_imfp.Ek, ele_imfp.imfp_S2_avg, '--', 'color', 'g', 'LineWidth', 2.0);
% - PLOTTING UNIVERSAL CURVE
loglog(ele_imfp.Ek, ele_imfp.imfp_universal, '-', 'color', 'k', 'LineWidth', 3.0);
% - FORMATTING THE FIGURE
% -- Plotting the x- and y-axes
a = line([0 0], [-1e5, 1e5], 'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
b = line([-1e5, 1e5], [0 0], 'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
a.Annotation.LegendInformation.IconDisplayStyle = 'off';
b.Annotation.LegendInformation.IconDisplayStyle = 'off';
% -- Labelling the x- and y-axes
xlabel('$$ \bf Electron\ Kinetic\ Energy\ [eV] $$', 'interpreter', 'latex');
ylabel('$$ \bf IMFP\ [Angstrom] $$', 'interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
ax.XLim = [1, 1e4]; ax.YLim = [1, 1e3];
ax.XScale = 'log'; ax.YScale = 'log';
text(0.05, 0.07, string(material), 'FontSize', 18, 'color', 'k', 'Units','normalized', 'FontWeight', 'bold');
% -- Add a legend
h       = zeros(6, 1);
h(1)    = plot(NaN,NaN,'k-', 'LineWidth', 3.0);
h(2)    = plot(NaN,NaN,'rx-');
h(3)    = plot(NaN,NaN,'b:', 'LineWidth', 2.0);
h(4)    = plot(NaN,NaN,'b--', 'LineWidth', 2.0);
h(5)    = plot(NaN,NaN,'g:', 'LineWidth', 2.0);
h(6)    = plot(NaN,NaN,'g--', 'LineWidth', 2.0);
legend(h, 'Universal','Optical/Experiment','TPP-2M','TPP-2M(avg)','S1','S2(avg)', 'fontsize', 10);
end