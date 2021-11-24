function bzone_props(bzOverlay, bz_args)
% bzone_props(bzOverlay, bz_args)
%   Function that will plot the BZ overlay once it has been determined via
%   the 'extract_lattice()' and 'extract_bz_overlay()' functions. Allows
%   the user to choose some properties of the BZ plots too. This function
%   allows for a constant, global formatting for Brilluoin Zone & 
%   Overlay properties within the PESTools package.
%
%   REQ. FUNCTIONS: none
%
%   IN:
%   -   bzOverlay:    	data structure of the sliced Iso-E ARPES data.
%   -   bz_args:        1x3 cell of {bz_offset [xShift,yShift], bz_col [1,1,1], bz_lwidth [1.5]}
%
%   OUT: (none) plot of the BZ overlay

%% Default parameters
if nargin < 2; bz_args = cell(1,3);  end
if isempty(bz_args); bz_args = cell(1,3); end

%% - 1 - Initialising the BZ overlay parameters
% - Extracting the BZ values
bz_offset   = bz_args{1}; if isempty(bz_offset); 	bz_offset = [0 0];  end
bz_col      = bz_args{2}; if isempty(bz_col);       bz_col = [1 1 1];   end
bz_lwidth   = bz_args{3}; if isempty(bz_lwidth);    bz_lwidth = 1.5;    end
% - Extracting the initial axis limits
xl = xlim; yl = ylim;

%% - 2 - Plotting the BZ overlay
for ii = 1:size(bzOverlay.X, 2)
    plot(bzOverlay.X{ii} - bz_offset(1), bzOverlay.Y{ii} - bz_offset(2), 'w-', 'color', bz_col, 'linewidth', bz_lwidth);
end
xticks(round(-20*bzOverlay.gX:bzOverlay.gX/2:20*bzOverlay.gX,2));
yticks(round(-20*bzOverlay.gY:bzOverlay.gY/2:20*bzOverlay.gY,2));
% - Pinning axes limits to original values
axis([xl(1), xl(2), yl(1), yl(2)]);
end