function hvOverlay = extract_hv_overlay(hv, navi_args, bzOverlay, plot_results)
% hvOverlay = extract_hv_overlay(hv, navi_args, bzOverlay, plot_results)
%   This is a function that extracts the crystallographic vectors
%   of a general crystal structure in both real- and reciprocal-
%   space based on the inputs. All 14 Bravais lattices can be 
%   defined, including primitive-, base-, face- and body-centered
%   lattices. The outputs are two MATLAB structures that contain
%   all of the real- and reciprocal-space data.
%
%   IN:
%   - hv:        	single or 1x2 vector of [minHv, maxHv] for the user defined photon energy.
%   - navi_args:	MATLAB data structure containing all of the navigation arguments below;
%       .(ePhi):                	Work function of analyser
%       .(alpha):                   Angle between analyser and sample normal
%       .(thtM):                	Primary manipulator rotation
%       .(eBref):                	Binding energy reference
%       .(thtAref):             	Analyser angle to sample normal in geometry is at origin
%       .(kxref):                  	kx reference is at origin
%       .(tltM):                    Primary tilt rotation is set to origin
%   	.(V000):                    User defined inner potential.
%   - bzOverlay:     MATLAB data structure containing all BZ plane information from 'extract_bz_overlay()';
%       .(crystal):             	string of the crystal type.
%       .(crystal_plane):        	string of the (h,k,l) plane extracted.
%       .(area):                 	area of the planar BZ unit cell.
%       .(X):                       cell array of the x vertices of the planar BZ cell slice.
%       .(Y):                       cell array of the y vertices of the planar BZ cell slice.
%       .(gX):                      double of the reciprocal-tesselation vector in x.
%       .(gY):                      double of the reciprocal-tesselation vector in y.
%   -  plot_results:         if 1, will plot figure of the results, otherwise it wont.

%% Default parameters
navi_default_args.ePhi      = 4.5; 
navi_default_args.alpha     = 20.0; 
navi_default_args.thtM      = 0; 
navi_default_args.eBref     = 0; 
navi_default_args.thtAref	= 0; 
navi_default_args.kxref     = 0; 
navi_default_args.tltM      = 0; 
navi_default_args.V000      = 12.57; 
if nargin < 2; navi_args = navi_default_args; bzOverlay = []; plot_results = 1; end
if nargin < 3; bzOverlay = []; plot_results = 1; end
if nargin < 4; plot_results = 1; end
if isempty(hv); hv = [350, 1500]; end
if isempty(navi_args); navi_args = navi_default_args; end
if isempty(plot_results); plot_results = 1; end
disp('-> Extracting k-space ARPES scan line...')

%% 1.1 - Extracting the wave-vector path in k-space for the ARPES scan
tht = linspace(-12, 12, 5e2);
if length(hv) == 1
    surfNormX = SurfNormX(hv, navi_args.eBref, navi_args.kxref, navi_args.thtM, navi_args.thtAref);
    kx = Kxx(hv, navi_args.eBref, navi_args.thtM, tht, surfNormX);
    kz = Kzz(hv, navi_args.eBref, navi_args.thtM, tht, navi_args.tltM, navi_args.V000, surfNormX);
elseif length(hv) == 2
    hv = sort(hv); hv = [hv(1), hv(end)];
    surfNormX(1,:) = SurfNormX(hv(1), navi_args.eBref, navi_args.kxref, navi_args.thtM, navi_args.thtAref);
    kx(1,:) = Kxx(hv(1), navi_args.eBref, navi_args.thtM, tht, surfNormX(1,:));
    kz(1,:) = Kzz(hv(1), navi_args.eBref, navi_args.thtM, tht, navi_args.tltM, navi_args.V000, surfNormX(1,:));
    surfNormX(2,:) = SurfNormX(hv(2), navi_args.eBref, navi_args.kxref, navi_args.thtM, navi_args.thtAref);
    kx(2,:) = Kxx(hv(2), navi_args.eBref, navi_args.thtM, tht, surfNormX(2,:));
    kz(2,:) = Kzz(hv(2), navi_args.eBref, navi_args.thtM, tht, navi_args.tltM, navi_args.V000, surfNormX(2,:));
end
%% 1.2 - Extracting the kz-broadening
eKE     = hv - navi_args.ePhi + navi_args.eBref;
imfp    = 10 .* (143 ./ eKE .^2 + 0.054.*sqrt(eKE));
dkz     = 1 ./ imfp;

%% 2 - Updating the figure with the photon energy path
if plot_results == 1
    figX = figure(360579); 
    figX.Name = 'BZ Navigation'; 
    hold on;
    % - Formatting the figure
    ax = gca; ax.TickLabelInterpreter = 'latex'; ax.TickDir = 'both';
    xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  k_{\perp} (\AA^{-1}) $$', 'Interpreter', 'latex');
    line([0, 0], [-1e5, 1e5], 'color', 'k', 'linestyle', '-');
    %% 2.1 - Plotting the BZ planar slice
    if ~isempty(bzOverlay)
        for i = 1:size(bzOverlay.X, 2)
            plot(bzOverlay.X{i}, bzOverlay.Y{i}, 'k-','linewidth', 2);
        end
        % -- Defining the axes properties
        xticks(round(-1e3*norm(bzOverlay.gX):0.5*norm(bzOverlay.gX):1e3*norm(bzOverlay.gX),2));
        yticks(round(-1e3*norm(bzOverlay.gY):0.5*norm(bzOverlay.gY):1e3*norm(bzOverlay.gY),2));
        title_txt = sprintf("%s - %s", bzOverlay.crystal, bzOverlay.crystal_plane);
        title(title_txt);
    end
    %% 2.2 - Plotting the k-space patches probed with ARPES
    col = rand(1,3);
    % - For a single photon energy
    if length(hv) == 1
        % -- Extract the planar limits of the ARPES scan
        kx_lower = (kx); kz_lower = (kz - dkz)';
        kx_upper = (kx); kz_upper = (kz + dkz)';
        % -- Plotting the ARPES scan plane which includes the Kz broadening
        patch([kx_upper, fliplr(kx_lower)], [kz_upper, fliplr(kz_lower)], col, 'edgecolor', 'none', 'facealpha', 0.5);
        % -- Plotting the ARPES scan line
        plot(kx, kz, '-', 'linewidth', 2, 'color', col);
        % -- Adding text
        text(max(kx(:)), min(kz_lower(:)), sprintf("%.1f eV", hv),...
            'color', col, 'Fontsize', 12, 'horizontalalignment', 'left');
    % - For a range of photon energies
    elseif length(hv) == 2
        % -- Extract the planar limits of the ARPES scan
        kx_lower = (kx(1,:)); kz_lower = (kz(1,:) - dkz(1));
        kx_upper = (kx(2,:)); kz_upper = (kz(2,:) + dkz(2));
        % -- Plotting the ARPES scan plane which includes the Kz broadening
        patch([kx_upper, fliplr(kx_lower)], [kz_upper, fliplr(kz_lower)], col, 'edgecolor', 'none', 'facealpha', 0.5);
        % -- Plotting the ARPES scan line
        plot(kx(1,:), kz(1,:), '-', 'linewidth', 2, 'color', col);
        plot(kx(2,:), kz(2,:), '-', 'linewidth', 2, 'color', col);
        % -- Adding text
        text(max(kx(1,:)), min(kz_lower(:)), sprintf("%.1f eV", hv(1)),...
            'color', col, 'Fontsize', 12, 'horizontalalignment', 'left');
        text(max(kx(2,:)), max(kz_upper(:)), sprintf("%.1f eV", hv(2)),...
            'color', col, 'Fontsize', 12, 'horizontalalignment', 'left');
    end
    % - Adding the axes limits
    axis equal; grid on;
    xlim([min(kx_upper(:)), max(kx_upper(:))]*1.05);
    ylim([0.95*min(kz_lower(:)), 1.05*max(kz_upper(:))]);
end
%% 4.0 - Assigning the data structure
hvOverlay               = struct();
hvOverlay.navi_args     = navi_args;
hvOverlay.bzOverlay     = bzOverlay;
hvOverlay.imfp          = imfp;
hvOverlay.dkz           = dkz;
hvOverlay.kx            = kx;
hvOverlay.kz            = kz;
% hvOverlay.kx_lower      = kx_lower;
% hvOverlay.kx_upper      = kx_upper;
end