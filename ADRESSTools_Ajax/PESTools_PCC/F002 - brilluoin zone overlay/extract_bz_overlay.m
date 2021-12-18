function bzOverlay = extract_bz_overlay(reciStr, crystal_plane, plot_results)
% bzOverlay = extract_bz_overlay(reciStr, crystal_plane, plot_results)
%   This is a function that extracts the planar Brilluoin zone 
%   along the defined crystal plane (h,k,l) in miller indices notation.
%   This Brilluoin zone plane then represents the planar slice that
%   that is probed through via ARPES by changing the scan parameter.
%
%   IN:
%   -  reciStr:         MATLAB data structure containing reciprocal-space data from 'extract_lattice()'
%   -  crystal_plane:   string of the crystal plane to extract in "(hkl)" format.
%   -  plot_results: 	if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   -  bzOverlay          MATLAB data structure containing all BZ plane information below;
%   	.(crystal):           string of the crystal type.
%   	.(crystal_plane):     string of the (h,k,l) plane extracted.
%		.(area):              area of the planar BZ unit cell.
%   	.(X):                 cell array of the x vertices of the planar BZ cell slice.
%   	.(Y):                 cell array of the y vertices of the planar BZ cell slice.
%		.(gX):                double of the reciprocal-tesselation vector in x.
%		.(gY):                double of the reciprocal-tesselation vector in y.

%% Default parameters
if nargin < 2; crystal_plane = "001"; plot_results = 0; end
if nargin < 3; plot_results = 0; end
if isempty(crystal_plane); crystal_plane = "(001)"; end
if isempty(plot_results); plot_results = 0; end
disp('-> Extracting plane slice of BZ...')
wbar = waitbar(0, 'Extracting BZ overlay...', 'Name', 'extract_bz_overlay');

%% Initialising the input variables
crystal = reciStr.crystal;
G1 =  reciStr.G1;
G2 =  reciStr.G2;
G3 =  reciStr.G3;
Gbz = reciStr.Gbz;
%% 1 - Extracting the First Brilluoin Zone boundary
% - Extracting the point cloud for the First Brilluoin Zone
x_pts = cell2mat(Gbz(:,1));
y_pts = cell2mat(Gbz(:,2));
z_pts = cell2mat(Gbz(:,3));
cloud_pts = [x_pts, y_pts, z_pts];
% - Projecting the point cloud onto the desired planar slice
% --- For the (100) plane
if crystal_plane == "(100)"
    proj_vec = [0,1,1];
    proj_pts = cloud_pts .* proj_vec;
    xpts1 = proj_pts(:,2);
    ypts1 = proj_pts(:,3);
% --- For the (010) plane
elseif crystal_plane == "(010)"
    proj_vec = [1,0,1];
    proj_pts = cloud_pts .* proj_vec;
    xpts1 = proj_pts(:,1);
    ypts1 = proj_pts(:,3);
% --- For the (001) plane
elseif crystal_plane == "(001)"
    proj_vec = [1,1,0];
    proj_pts = cloud_pts .* proj_vec;
    xpts1 = proj_pts(:,1);
    ypts1 = proj_pts(:,2);
% --- For the (110) plane
elseif crystal_plane == "(110)"
    proj_vec = [1,0,1];
    R = [cosd(45) -sind(45) 0; sind(45) cosd(45) 0; 0 0 1];
    cloud_pts = cloud_pts * R;
    proj_pts = cloud_pts .* proj_vec;
    xpts1 = proj_pts(:,1);
    ypts1 = proj_pts(:,3);
end
% - Shifting the boundary points to the origin
xpts1 = xpts1 - (max(xpts1(:))+min(xpts1(:)))/2;
ypts1 = ypts1 - (max(ypts1(:))+min(ypts1(:)))/2;

%% 2 - Creating the boundary Brilluoin Zone polygon and finding the area of a single zone
[b_indx, area] = boundary(xpts1, ypts1, 0.1);
xpts1 = xpts1(b_indx);
ypts1 = ypts1(b_indx);

%% 3 - Tessalating the polygons to get the overlay
l=1;
% Iterating over all (or most of!) reciprocal space
for i = -25:25
    for j = -25:25
        % -- For the (100) plane
        if crystal_plane == "(100)"
            % -- Applying the tessalated shifts
            X{l} = xpts1 + i*G2(2) + j*G3(2);
            Y{l} = ypts1 + i*G2(3) + j*G3(3);
        % -- For the (010) plane
        elseif crystal_plane == "(010)"
            % -- Applying the tessalated shifts
            X{l} = xpts1 + i*G1(1) + j*G3(1);
            Y{l} = ypts1 + i*G1(3) + j*G3(3);
        % -- For the (001) plane
        elseif crystal_plane == "(001)"
            % -- Applying the tessalated shifts
            X{l} = xpts1 + i*G1(1) + j*G2(1);
            Y{l} = ypts1 + i*G1(2) + j*G2(2);
        % -- For the (110) plane
        elseif crystal_plane == "(110)"
            % --- Extracting the shift vectors
            xshift = norm(0.5*G1 + 0.5*G2);
            yshift = norm(0.5*G3);
            % -- Applying the tessalated shifts
            if mod(i,2) == 1
                X{l} = xpts1 + i*xshift;
                Y{l} = ypts1 + yshift +  j*2*max(yshift);
            else
                X{l} = xpts1 + i*xshift;
                Y{l} = ypts1 + j*2*max(yshift);
            end
        end
        l = l + 1;  
    end
end
%% Assigning the overlay structure to a MATLAB structure
waitbar(1, wbar, 'Assigning the BZ overlay variables...', 'Name', 'extract_lattice');
% - Real-space information
bzOverlay.crystal = crystal;
bzOverlay.crystal_plane = crystal_plane;
bzOverlay.area = area;
bzOverlay.X = X;
bzOverlay.Y = Y;
if crystal_plane == "(100)"
    bzOverlay.gX = norm(G2);
    bzOverlay.gY = norm(G3);
elseif crystal_plane == "(010)"
    bzOverlay.gX = norm(G1);
    bzOverlay.gY = norm(G3);
elseif crystal_plane == "(001)"
    bzOverlay.gX = norm(G1);
    bzOverlay.gY = norm(G2);
elseif crystal_plane == "(110)"
    bzOverlay.gX = norm(G1 + G2);
    bzOverlay.gY = norm(G3);
end
%% Close wait-bar
close(wbar);
%% Figure summary of the real- and reciprocal-space structures
if plot_results == 1

    fig = figure(); set(fig, 'position', [1,1,950,650]);

    %% - RECIPROCAL-SPACE FIGURE
    subplot(1,2,1); hold on;
    % -- Plotting the axes lines
    line([min(reciStr.G(:,1)), max(reciStr.G(:,1))], [0 0], [0,0],  'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
    line([0 0], [min(reciStr.G(:,2)), max(reciStr.G(:,2))], [0, 0], 'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
    line([0 0], [0, 0], [min(reciStr.G(:,3)), max(reciStr.G(:,3))], 'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
    % -- Plotting the translational real-space lattice vectors
    plot3([0,reciStr.G1(1)],[0,reciStr.G1(2)],[0,reciStr.G1(3)], '-', 'linewidth', 3, 'color', [0.3, 0.8, 0.3]);
    plot3([0,reciStr.G2(1)],[0,reciStr.G2(2)],[0,reciStr.G2(3)], '-', 'linewidth', 3, 'color', [0.3, 0.8, 0.3]);
    plot3([0,reciStr.G3(1)],[0,reciStr.G3(2)],[0,reciStr.G3(3)], '-', 'linewidth', 3, 'color', [0.3, 0.8, 0.3]);
    % -- Plotting the multiple translations of the real-space vectors
    plot3(reciStr.G(:,1),reciStr.G(:,2),reciStr.G(:,3),'b.','markersize',20);
    % -- Plotting all the nearest neighbour points
    for i = 1:2:size(reciStr.Gnn,1)
        plot3([reciStr.Gnn(i,1), reciStr.Gnn(i+1,1)], [reciStr.Gnn(i,2), reciStr.Gnn(i+1,2)], [reciStr.Gnn(i,3), reciStr.Gnn(i+1,3)], 'k-', 'linewidth', 1);
    end
    % -- Plotting the first Brilluoin Zone
    for i = 1:size(reciStr.Gbz,1)
        patch(reciStr.Gbz{i,1}, reciStr.Gbz{i,2}, reciStr.Gbz{i,3},[0.3 0.8 0.3], 'FaceAlpha',0.75, 'EdgeColor', 'none');
    end
    % -- Plotting a path of the brilluoin zone slice
    xP = [-1e3, 1e3, 1e3, -1e3, -1e3];
    yP = [1e3, 1e3, -1e3, -1e3, 1e3];
    zP = [0, 0, 0, 0, 0];
    if crystal_plane == "(100)"
        patch(zP, xP, yP, [0.7, 0.7, 0.7], 'FaceAlpha',0.75, 'EdgeColor', 'none');
    elseif crystal_plane == "(010)"
        patch(xP, zP+reciStr.G2(2), yP, [0.7, 0.7, 0.7], 'FaceAlpha',0.75, 'EdgeColor', 'none');
    elseif crystal_plane == "(001)"
        patch(xP, yP, zP, [0.7, 0.7, 0.7], 'FaceAlpha',0.75, 'EdgeColor', 'none');
    elseif crystal_plane == "(110)"
        xP = [1e3, 1e3, -1e3, -1e3, 1e3];
        yP = [1e3, 1e3, -1e3, -1e3, 1e3]+reciStr.G2(2);
        zP = [-1e3, 1e3, 1e3, -1e3, -1e3];
        patch(xP, yP, zP, [0.7, 0.7, 0.7], 'FaceAlpha',0.75, 'EdgeColor', 'none');
    end
    % -- Defining the axes properties
    % - Figure formatting
    ax = gca;
    % Font properties
    ax.FontName = 'Helvetica'; ax.FontSize = 14;
    % Tick properties
    ax.TickLabelInterpreter = 'latex';
    ax.TickDir = 'both';
    % Box Styling properties
    ax.LineWidth = 1.2;
    % Axis labels, limits and ticks
    xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex');
    zlabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex');
    xticks(round(-10*norm(G1):norm(G1):10*norm(G1),2));
    yticks(round(-10*norm(G2):norm(G2):10*norm(G2),2));
    zticks(round(-10*norm(G3):norm(G3):10*norm(G3),2));
    % -- Modify the view
    view(3); camlight(-30,24);
    title_txt2 = sprintf("%s; Reciprocal; Brilluoin Zone", crystal);
    title(title_txt2);
    axis tight equal; rotate3d on;
    pbaspect([1,1,1]);
    axis(reciStr.Glims);

    %% - BRILLUOIN ZONE PLANAR SLICE
    subplot(1,2,2); hold on;
    for i = 1:size(X, 2)
        plot(X{i}, Y{i}, 'k-','linewidth', 2);
    end
    % -- Defining the axes properties
    % - Figure formatting
    ax = gca;
    % Font properties
    ax.FontName = 'Helvetica'; ax.FontSize = 14;
    % Tick properties
    ax.TickLabelInterpreter = 'latex';
    ax.TickDir = 'both';
    % Box Styling properties
    ax.LineWidth = 1.2;
    % Axis labels, limits and ticks
    title_txt = sprintf("%s; BZNavi through %s", bzOverlay.crystal, bzOverlay.crystal_plane);
    title(title_txt);
    xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex');
    xticks(round(-100*norm(bzOverlay.gX):0.5*norm(bzOverlay.gX):100*norm(bzOverlay.gX),2));
    yticks(round(-100*norm(bzOverlay.gY):0.5*norm(bzOverlay.gY):100*norm(bzOverlay.gY),2));
    axis(1.5*[-norm(bzOverlay.gX), norm(bzOverlay.gX), -norm(bzOverlay.gY),norm(bzOverlay.gY)]);
    axis equal;
    % Axis labels, limits and ticks
    title_txt = sprintf("BZ plane; %s", crystal_plane);
    title(title_txt);
    axis equal;
    if crystal_plane == "(100)"
        xlabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex');
        ylabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex');
    elseif crystal_plane == "(010)"
        xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
        ylabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex');
    elseif crystal_plane == "(001)"
        xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
        ylabel('$$ \bf  k_y (\AA^{-1}) $$', 'Interpreter', 'latex');
    elseif crystal_plane == "(110)"
        xlabel('$$ \bf  k_x (\AA^{-1}) $$', 'Interpreter', 'latex');
        ylabel('$$ \bf  k_z (\AA^{-1}) $$', 'Interpreter', 'latex');
    end
    xticks(round(-10*bzOverlay.gX:0.5*bzOverlay.gX:10*bzOverlay.gX,2));
    yticks(round(-10*bzOverlay.gY:0.5*bzOverlay.gY:10*bzOverlay.gY,2));
    axis([-bzOverlay.gX, bzOverlay.gX, -bzOverlay.gY, bzOverlay.gY]*1.05);
end
end