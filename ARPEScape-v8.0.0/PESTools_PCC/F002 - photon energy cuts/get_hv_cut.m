function hvCut = get_hv_cut(hv, navi_args, plot_results)
% hvCut = get_hv_cut(hv, navi_args, plot_results)
%   This is a function that extracts the ARPES cut made in reciprocal space
%   for a defined measurement geometry at a given photon energy. 
%
%   IN:
%   - hv:        	single or 1x2 vector of [minHv, maxHv] for the user defined photon energy.
%   - navi_args:	MATLAB data structure containing all of the navigation arguments below;
%       .(ePhi):                	Work function of analyser
%       .(alpha):                   Angle between analyser and sample normal (July 2022: 20deg, from August 2022: 9deg).
%       .(thtM):                	Primary manipulator rotation
%       .(eBref):                	Binding energy reference
%       .(thtAref):             	Analyser angle to sample normal in geometry is at origin
%       .(kxref):                  	kx reference is at origin
%       .(tltM):                    Primary tilt rotation is set to origin
%   	.(V000):                    User defined inner potential.
%   -  plot_results: 	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -  hvCut        MATLAB data structure containing all hv cut information below;
%   	.(hv):          input photon energy
%   	.(navi_args):   input navigation arguments
%		.(kx):          500x1 (or 500x2) matrix of the cut kx values.
%   	.(kz):          500x1 (or 500x2) matrix of the cut kz values.

%% Default parameters
navi_default_args.ePhi      = 4.5; 
navi_default_args.alpha     = 9.0;  % present geometry 9deg; Old geometry 20deg;
navi_default_args.thtM      = 0; 
navi_default_args.eBref     = 0; 
navi_default_args.thtAref	= 0; 
navi_default_args.kxref     = 0; 
navi_default_args.tltM      = 0; 
navi_default_args.V000      = 12.57; 
if nargin < 2; navi_args = navi_default_args; plot_results = 0; end
if nargin < 3; plot_results = 0; end
if isempty(navi_args); navi_args = navi_default_args; end
if isempty(plot_results); plot_results = 0; end

%% Initialising the input variables
tht = linspace(-12, 12, 5e2);
%% 1 - Extracting the wave-vector path in k-space for the ARPES scan
if length(hv) == 1
    surfNormX(:,1)      = SurfNormX(hv, navi_args.eBref, navi_args.kxref, navi_args.thtM, navi_args.thtAref, navi_args.alpha);
    kx(:,1)             = Kxx(hv, navi_args.eBref, navi_args.thtM, tht, surfNormX, navi_args.alpha);
    kz(:,1)             = Kzz(hv, navi_args.eBref, navi_args.thtM, tht, navi_args.tltM, navi_args.V000, surfNormX, navi_args.alpha);
elseif length(hv) == 2
    hv = sort(hv); hv = [hv(1), hv(end)];
    % -- first photon energy
    surfNormX(1,:)      = SurfNormX(hv(1), navi_args.eBref, navi_args.kxref, navi_args.thtM, navi_args.thtAref, navi_args.alpha);
    kx(:,1)             = Kxx(hv(1), navi_args.eBref, navi_args.thtM, tht, surfNormX(1,:), navi_args.alpha);
    kz(:,1)             = Kzz(hv(1), navi_args.eBref, navi_args.thtM, tht, navi_args.tltM, navi_args.V000, surfNormX(1,:), navi_args.alpha);
    % -- second photon energy
    surfNormX(2,:)      = SurfNormX(hv(2), navi_args.eBref, navi_args.kxref, navi_args.thtM, navi_args.thtAref, navi_args.alpha);
    kx(:,2)             = Kxx(hv(2), navi_args.eBref, navi_args.thtM, tht, surfNormX(2,:), navi_args.alpha);
    kz(:,2)             = Kzz(hv(2), navi_args.eBref, navi_args.thtM, tht, navi_args.tltM, navi_args.V000, surfNormX(2,:), navi_args.alpha);
elseif length(hv) > 2
    error('hv can only be a double single or 1x2 vector of [minHv, maxHv].');
end
%% 2 - Assigning the hv cut to MATLAB data structure
hvCut               = struct();
hvCut.hv            = hv;
hvCut.navi_args     = navi_args;
hvCut.kx            = kx;
hvCut.kz            = kz;
%% 3 - Figure summary of the hv cuts
if plot_results == 1
    fig = figure(); set(fig, 'position', [1,1,600,600]);
    hold on;
    col = lines(2);
    % - For a single photon energy
    if length(hv) == 1
        % -- Plotting the ARPES scan line
        plot(kx, kz, '-', 'linewidth', 2, 'color', col(1,:));
        % -- Adding text
        text(max(kx(:)), min(kz(:)), sprintf("%.0f eV", hv),...
            'color', col(1,:), 'Fontsize', 12, 'horizontalalignment', 'left', 'FontWeight','bold');
    % - For a range of photon energies
    elseif length(hv) == 2
        % -- Plotting the ARPES scan plane which includes the Kz broadening
        patch([kx(:,2); flipud(kx(:,1))], [kz(:,2); flipud(kz(:,1))], col(2,:), 'edgecolor', 'none', 'facealpha', 0.5);
        % -- Plotting the ARPES scan line
        plot(kx(:,1), kz(:,1), '-', 'linewidth', 2, 'color', col(1,:));
        plot(kx(:,2), kz(:,2), '-', 'linewidth', 2, 'color', col(2,:));
        % -- Adding text
        text(max(kx(:,1)), min(kz(:,1)), sprintf("%.0f eV", hv(1)),...
            'color', col(1,:), 'Fontsize', 12, 'horizontalalignment', 'left', 'FontWeight','bold');
        text(max(kx(:,2)), min(kz(:,2)), sprintf("%.0f eV", hv(2)),...
            'color', col(2,:), 'Fontsize', 12, 'horizontalalignment', 'left', 'FontWeight','bold');
    end
    % - Adding the axes limits
    axis equal; grid on;
    xlim([min(kx(:)), max(kx(:))]*1.05);
    ylim([0.95*min(kz(:)), 1.05*max(kz(:))]);
end
end