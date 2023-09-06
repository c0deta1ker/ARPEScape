function hvCut = get_hv_cut(hv, EPhi, V000, ThtM, TltM, InciAlpha, RefEB, RefThtA, RefKx, plot_results)
% hvCut = get_hv_cut(hv, EPhi, V000, ThtM, TltM, InciAlpha, RefEB, RefThtA, RefKx, plot_results)
%   This is a function that extracts the ARPES cut made in reciprocal space
%   for a defined measurement geometry at a given photon energy. 
%
%   IN:
%   - hv:        	User defined photon energy (eV)
%   - EPhi:         Work function of analyser (eV) (@ADRESS: 4.50 eV)
%   - V000:         User defined inner potential (eV)
%   - ThtM:         Primary manipulator rotation (degree) (at normal emission = 0)
%   - TltM:         Primary tilt rotation (degree) (at normal emission = 0)
%   - InciAlpha:    Angle between analyser and sample normal (@ADRESS: < July 2022 = 20deg, present = 9deg).
%   - RefEB:        Reference Binding Energy (eV)
%   - RefThtA:      Reference Analyser Angle relative to sample normal (degree) (at normal emission = 0)
%   - RefKx:        Reference Kx value at the angle Ref_ThtA (Angstrom^-1)
%   - plot_results: if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -  hvCut        MATLAB data structure containing all hv cut information below;
%   	.(hv):          input photon energy
%   	.(navi_args):   input navigation arguments
%		.(kx):          500x1 (or 500x2) matrix of the cut kx values.
%   	.(kz):          500x1 (or 500x2) matrix of the cut kz values.

%% Default parameters
if nargin < 2; EPhi = 4.50; end
if nargin < 3; V000 = 12.57; end
if nargin < 4; ThtM = 0.00; end
if nargin < 5; TltM = 0.00; end
if nargin < 6; InciAlpha = 9.00; end
if nargin < 7; RefEB = 0.00; end
if nargin < 8; RefThtA = 0.00; end
if nargin < 9; RefKx = 0.00; end
if nargin < 10; plot_results = 0; end
if isempty(EPhi);       EPhi = 4.50; end
if isempty(V000);       V000 = 12.57; end
if isempty(ThtM);       ThtM = 0.00; end
if isempty(TltM);       TltM = 0.00; end
if isempty(InciAlpha);  InciAlpha = 9.00; end
if isempty(RefEB);      RefEB = 0.00; end
if isempty(RefThtA);    RefThtA = 0.00; end
if isempty(RefKx);      RefKx = 0.00; end
if isempty(plot_results); plot_results = 0; end

%% Initialising the input variables
tht = linspace(-12, 12, 5e2);
%% 1 - Extracting the wave-vector path in k-space for the ARPES scan
if length(hv) == 1
    surfNormX(:,1)      = SurfNormX(hv, RefEB, RefKx, ThtM, RefThtA, InciAlpha);
    kx(:,1)             = Kxx(hv, RefEB, ThtM, tht, surfNormX, InciAlpha);
    kz(:,1)             = Kzz(hv, RefEB, ThtM, tht, TltM, V000, surfNormX, InciAlpha);
else
    error('hv can only be a single scalar value.');
end
%% 2 - Assigning the hv cut to MATLAB data structure
hvCut                   = struct();
hvCut.hv                = hv;
hvCut.args.EPhi         = EPhi;
hvCut.args.V000         = V000;
hvCut.args.ThtM         = ThtM;
hvCut.args.TltM         = TltM;
hvCut.args.InciAlpha    = InciAlpha;
hvCut.args.RefEB        = RefEB;
hvCut.args.RefThtA      = RefThtA;
hvCut.args.RefKx        = RefKx;
hvCut.kx                = kx;
hvCut.kz                = kz;
%% 3 - Figure summary of the hv cuts
if plot_results == 1
    fig = figure(); set(fig, 'position', [1,1,600,600]);
    hold on;
    % -- Plotting the ARPES scan line
    plot(kx, kz, '-', 'linewidth', 2, 'color', 'b');
    % -- Adding text
    text(max(kx(:)), min(kz(:)), sprintf("%.0f eV", hv), 'color', 'b', 'Fontsize', 12, 'horizontalalignment', 'left', 'FontWeight','bold');
    % - Adding the axes limits
    axis equal; grid on;
    xlim([min(kx(:)), max(kx(:))]*1.05);
    ylim([0.95*min(kz(:)), 1.05*max(kz(:))]);
end
end