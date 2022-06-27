function [realStr, reciStr] = extract_lattice(crystal_args, plot_results)
% [realStr, reciStr] = extract_lattice(crystal_args, plot_results)
%   This is a function that extracts the crystallographic vectors
%   of a general crystal structure in both real- and reciprocal-
%   space based on the inputs. All 14 Bravais lattices can be 
%   defined, including primitive-, base-, face- and body-centered
%   lattices. The outputs are two MATLAB structures that contain
%   all of the real- and reciprocal-space data.
%
%   IN:
%   - crystal_args:   1x3 cell of {crystal, (a, b, c), (alpha, beta, gamma)} described below;
%           -- crystalType;             "CUB-cP-Oh",  "CUB-cI-Oh", "CUB-cF-Oh", 
%                                       "HEX-hP-D6h", "RHL-hR-D3d", 
%                                       "TET-tP-D4h", "TET-tI-D4h", 
%                                       "ORC-oP-D2h", "ORC-oS-D2h", "ORC-oI-D2h", "ORC-oF-D2h", 
%                                       "MCL-mP-C2h", "MCL-mS-C2h",
%                                       "TRI-aP-Ci" .
%           -- (a, b, c);            1x3 row-vector of the side lengths of the crystal systems unit cell (Angstroms).
%           -- (alpha, beta, gamma); 1x3 row vector the opposite angles of the crystal systems unit cell (degrees).
%   - plot_results:         if 1, will plot figure of the results, otherwise it wont.
%
%   OUT:
%   - realStr: MATLAB data structure containing all real-space information below.
%   	.(crystal):     string of the crystal type.
%   	.(T1):          a1 primitive translation vector.
%		.(T2):          a2 primitive translation vector.
%   	.(T3):          a3 primitive translation vector.
%   	.(Tb0):         atomic basis vectors.
%		.(Tvol):        volume of the real-space unit cell.
%		.(T):           atomic positions over multiple primitive-translations to show crystal structure.
%		.(Tlims):       axes limits in 3D for the lattice structure T.
%		.(Tnn):         linked lines between nearest neighbour atomsdisplayed in .(T).
%		.(Tws):         Wigner-Seitz cell in real-space.
%   - reciStr: MATLAB data structure containing all reciprocal-space information below.
%   	.(crystal):     string of the crystal type.
%   	.(G1):          b1 primitive reciprocal vector.
%		.(G2):          b2 primitive reciprocal vector.
%   	.(G3):          b3 primitive reciprocal vector.
%   	.(Gb0):         atomic basis vectors in reciprocal space.
%		.(Gvol):        volume of the reciprocal-space unit cell.
%		.(G):           reciprocal lattice over multiple translations.
%		.(Glims):       axes limits in 3D for the lattice structure G.
%		.(Gnn):         linked lines between nearest neighbour reciprocal-points displayed in .(G).
%		.(Gbz):         First Brilluoin zone in reciprocal-space.

%% Default parameters
if nargin < 1; crystal_args = {"CUB-cF-Oh", [5.431, 5.431, 5.431], [90, 90, 90]}; plot_results = 0; end
if nargin < 2; plot_results = 0; end
if isempty(crystal_args); crystal_args = {"CUB-cF-Oh", [5.431, 5.431, 5.431], [90, 90, 90]};  end
if isempty(plot_results); plot_results = 0;  end
disp('-> Extracting 3D BZ...')
wbar = waitbar(0, 'Extracting crystallographic variables...', 'Name', 'extract_lattice');

%% Initialising the input variables
crystal     = crystal_args{1};
a           = crystal_args{2}(1);
b           = crystal_args{2}(2);
c           = crystal_args{2}(3);
alpha       = crystal_args{3}(1) * pi/180;
beta        = crystal_args{3}(2) * pi/180;
gamma       = crystal_args{3}(3) * pi/180;

%% 1 - Defining the real-space unit cell from a general triclinic geometry
% -- Extracting the triclinic coefficients for real-space unit vectors
w1 = cos(alpha) - cos(beta) * cos(gamma)/(sin(beta) * sin(gamma));
w2 = sin(gamma)^2 - cos(beta)^2 - cos(alpha)^2 + 2*cos(alpha) *cos(beta) * cos(gamma);
w2 = sqrt(w2) / (sin(beta) * sin(gamma));
% -- Defining the real-space unit vectors
T1 = [a, 0, 0];
T2 = [b * cos(gamma), b	* sin(gamma), 0];
T3 = [c * cos(beta), c * w1*sin(beta), c*w2*sin(beta)];
% -- Extracting the volume of the real-space unit cell
Tvol = cross(T1, T2)*T3';
%% 2 - Defining the basis of the Bravais Lattice out of the 14 possible choices
% - 2.1 - The 7 Primitive (Simple) Bravais lattices
if crystal == "CUB-cP-Oh" || crystal == "HEX-hP-D6h" || crystal == "RHL-hR-D3d" ||...
        crystal == "TET-tP-D4h" || crystal == "ORC-oP-D2h" || crystal == "MCL-mP-C2h" ||...
        crystal == "TRI-aP-Ci"
    Tb0 = [0, 0, 0];
% - 2.2 - The 2 Base-Centered Bravais Lattices
elseif crystal == "ORC-oS-D2h" || crystal == "MCL-mS-C2h"
    Tb0 = [0, 0, 0;...
        0.5*(T1+T2)];
% - 2.3 - The 3 Body-Centered Bravais Lattices
elseif crystal == "CUB-cI-Oh" || crystal == "TET-tI-D4h" || crystal == "ORC-oI-D2h"
    Tb0 = [0, 0, 0;...
        0.5*(T1+T2+T3)];
% - 2.4 - The 2 Face-Centered Bravais Lattices
elseif crystal == "CUB-cF-Oh" || crystal == "ORC-oF-D2h"
    Tb0 = [0, 0, 0;...
        0.5*(T1+T2);...
        0.5*(T1+T3);...
        0.5*(T2+T3)];
end
%% 3 - Extracting the real-space lattice over multiple translations of T
l=1;
for i=-1:1
    for j=0:2
        for k=-1:1
            for n = 1:size(Tb0,1)
                T(l,:)=i*T1+j*T2+k*T3+Tb0(n,:);
                l=l+1;
            end
        end
    end
end
% - 3.1 - Finding the axes limited for a 2x2 super structure as in T
if crystal == "CUB-cP-Oh" || crystal == "HEX-hP-D6h" || crystal == "RHL-hR-D3d" ||...
        crystal == "TET-tP-D4h" || crystal == "ORC-oP-D2h" || crystal == "MCL-mP-C2h" ||...
        crystal == "TRI-aP-Ci"
    Tlims = 1.04*[min(T(:,1)), max(T(:,1)), min(T(:,2)), max(T(:,2)), min(T(:,3)), max(T(:,3))];
elseif crystal == "ORC-oS-D2h" || crystal == "MCL-mS-C2h" || crystal == "CUB-cI-Oh" ||...
        crystal == "TET-tI-D4h" || crystal == "ORC-oI-D2h" || crystal == "CUB-cF-Oh" ||...
        crystal == "ORC-oF-D2h" 
    Tlims = 1.04*[-1*norm(T1), norm(T1), 0, 2*norm(T2),-1*norm(T3), norm(T3)];
end
%% 4 - Finding all vector paths to link all nearest neighbour points
Tnn = [];
% - Filing through all possible points given in the T-matrix
for i = 1:length(T)-1
    waitbar(i/length(T)-1, wbar, 'Extracting real-space lattice...', 'Name', 'extract_lattice');
    for j = i+1:length(T)
        % -- Finding the Euclidean distance between the two lattice points
        dist  = sqrt((T(i,1)-T(j,1)).^2+(T(i,2)-T(j,2)).^2+(T(i,3)-T(j,3)).^2);
        %% 4.1 - Applying geometrical constraints for nearest neighbour points

        % --- For the 7 Primitive (Simple) Bravais lattices
        % Points are linked if their distance is identical to the
        % length of T1, T2 or T3.
        if crystal_args{1} == "CUB-cP-Oh" || crystal_args{1} == "HEX-hP-D6h" || crystal_args{1} == "RHL-hR-D3d" ||...
                crystal_args{1} == "TET-tP-D4h" || crystal_args{1} == "ORC-oP-D2h" || crystal_args{1} == "MCL-mP-C2h" ||...
                crystal_args{1} == "TRI-aP-Ci"
            if (abs(dist-norm(T1))) <= 0.0001 || (abs(dist-norm(T2))) <= 0.0001 || (abs(dist-norm(T3))) <= 0.0001                
                Tnn = vertcat(Tnn, [[T(i,1);T(j,1)], [T(i,2);T(j,2)], [T(i,3);T(j,3)]]);
            end

        % --- For the 2 Base-Centered Bravais Lattices
        % Points are linked if their distance is equal to the length
        % of the BCC basis vector. Also, if their length is equal to 
        % T1, T2 or T3, with the condition that it is not divisible by
        % the BCC basis vector (as the base-centered points are only
        % joined within the horizontal plane of the vertices).
        elseif crystal_args{1} == "ORC-oS-D2h" || crystal_args{1} == "MCL-mS-C2h" 
            if abs(dist-norm(Tb0(2,:))) <= 0.0001
                Tnn = vertcat(Tnn, [[T(i,1);T(j,1)], [T(i,2);T(j,2)], [T(i,3);T(j,3)]]);
            elseif (mod(T(i,2), 2*Tb0(2,2)) == 0 ||  abs(T(i,2)) <= 0.0001 ) && ((abs(dist-norm(T1))) <= 0.0001 || (abs(dist-norm(T2))) <= 0.0001 || (abs(dist-norm(T3))) <= 0.0001)
                Tnn = vertcat(Tnn, [[T(i,1);T(j,1)], [T(i,2);T(j,2)], [T(i,3);T(j,3)]]);
            end

        % --- For the 3 Body-Centered Bravais Lattices
        % Points are linked if their distance is equal to the length
        % of the BCC basis vector. Also, if their length is equal to 
        % T1, T2 or T3, with the condition that it is not divisible by
        % the BCC basis vector (as the body-centered points are only
        % joined radially outwards to vertices of unit cell).
        elseif crystal_args{1} == "CUB-cI-Oh" || crystal_args{1} == "TET-tI-D4h" || crystal_args{1} == "ORC-oI-D2h"
            if abs(dist-norm(Tb0(2,:))) <= 0.0001
                Tnn = vertcat(Tnn, [[T(i,1);T(j,1)], [T(i,2);T(j,2)], [T(i,3);T(j,3)]]);
            elseif mod(T(i,3), 2*Tb0(2,3)) == 0 && ((abs(dist-norm(T1))) <= 0.0001 || (abs(dist-norm(T2))) <= 0.0001 || (abs(dist-norm(T3))) <= 0.0001)
                Tnn = vertcat(Tnn, [[T(i,1);T(j,1)], [T(i,2);T(j,2)], [T(i,3);T(j,3)]]);
            end

        % --- For the 2 Face-Centered Bravais Lattices
        % Points are linked if their distance is equal to the length
        % of the FCC basis vector, with the proviso that the two points
        % are coplanar and not an FCC point. This is applied along each
        % dimension. Furthermore, if their length is equal to T1, T2 or T3, 
        % with the condition that it is not divisible by the FCC basis vector
        % (as the face-centered points are only joined radially outwards to
        % vertices of unit cell, not to other FCC points).
        elseif crystal_args{1} == "CUB-cF-Oh" || crystal_args{1} == "ORC-oF-D2h"
            if abs(dist-norm(Tb0(2,:))) <= 0.001 && abs(T(i,1)-T(j,1)) <= 0.001 && (abs(T(i,1)) <= 0.001 || mod(T(i,1), 2*Tb0(2,1)) == 0)
                Tnn = vertcat(Tnn, [[T(i,1);T(j,1)], [T(i,2);T(j,2)], [T(i,3);T(j,3)]]);
            elseif abs(dist-norm(Tb0(2,:))) <= 0.001 && abs(T(i,2)-T(j,2)) <= 0.001 && (abs(T(i,2)) <= 0.001 || mod(T(i,2), 2*Tb0(2,2)) == 0)
                Tnn = vertcat(Tnn, [[T(i,1);T(j,1)], [T(i,2);T(j,2)], [T(i,3);T(j,3)]]);
            elseif abs(dist-norm(Tb0(3,:))) <= 0.001 && abs(T(i,3)-T(j,3)) <= 0.001 && (abs(T(i,3)) <= 0.001 || mod(T(i,3), 2*Tb0(3,3)) == 0)
                Tnn = vertcat(Tnn, [[T(i,1);T(j,1)], [T(i,2);T(j,2)], [T(i,3);T(j,3)]]);
            elseif (abs(dist-norm(T1))) <= 0.0001 && abs(T(i,1)-T(j,1)) <= 0.001 && (abs(T(i,1)) <= 0.001 || mod(T(i,1), 2*Tb0(2,1)) == 0) && (abs(T(i,3)) <= 0.001 || mod(T(i,3), 2*Tb0(3,3)) == 0)
                Tnn = vertcat(Tnn, [[T(i,1);T(j,1)], [T(i,2);T(j,2)], [T(i,3);T(j,3)]]);
            elseif (abs(dist-norm(T2))) <= 0.0001 && abs(T(i,2)-T(j,2)) <= 0.001 && (abs(T(i,2)) <= 0.001 || mod(T(i,2), 2*Tb0(2,2)) == 0) && (abs(T(i,3)) <= 0.001 || mod(T(i,3), 2*Tb0(3,3)) == 0)
                Tnn = vertcat(Tnn, [[T(i,1);T(j,1)], [T(i,2);T(j,2)], [T(i,3);T(j,3)]]);
            end

        end
    end
end
%% 5 - Computing the Voronoi diagram, which is identical to the Wigner-Seitz cell
% - Calculating Voronoi diagram
[T_c, T_v] = voronoin(T); 
T_nx = T_c(T_v{floor(l/2)},:);
T_cell = convhulln(T_nx);
% - Defining the T_ws variable and shifting Wigner-Seitz cell to +T2
for i = 1:size(T_cell,1)
    Tws{i,1} = T_nx(T_cell(i,:),1)-T(floor(l/2),1)+T2(1);
    Tws{i,2} = T_nx(T_cell(i,:),2)-T(floor(l/2),2)+T2(2);
    Tws{i,3} = T_nx(T_cell(i,:),3)-T(floor(l/2),3)+T2(3);
end


%% A - Defining the reciprocal-space unit cell
% - A.1 - For simple and base-centered lattices, use T-vectors only
if crystal_args{1} == "CUB-cP-Oh" || crystal_args{1} == "HEX-hP-D6h" || crystal_args{1} == "RHL-hR-D3d" ||...
        crystal_args{1} == "TET-tP-D4h" || crystal_args{1} == "ORC-oP-D2h" || crystal_args{1} == "MCL-mP-C2h" ||...
        crystal_args{1} == "TRI-aP-Ci" || crystal_args{1} == "ORC-oS-D2h" || crystal_args{1} == "MCL-mS-C2h"
    G1 = 2*pi*cross(T2,T3)/Tvol;
    G2 = 2*pi*cross(T3,T1)/Tvol;
    G3 = 2*pi*cross(T1,T2)/Tvol;
% - A.2 - For body- and face-centered lattices, scale by a factor of 2
elseif crystal_args{1} == "CUB-cI-Oh" || crystal_args{1} == "TET-tI-D4h" || crystal_args{1} == "ORC-oI-D2h" ||...
        crystal_args{1} == "CUB-cF-Oh" || crystal_args{1} == "ORC-oF-D2h"
    G1 = 4*pi*cross(T2,T3)/Tvol;
    G2 = 4*pi*cross(T3,T1)/Tvol;
    G3 = 4*pi*cross(T1,T2)/Tvol;
end
% -- Extracting the volume of reciprocal unit cell
Gvol = cross(G1,G2)*G3';
%% B - Extracting the corresponding Brilluoin Zone in reciprocal-space
% - B.1 - The 7 Primitive (Simple) Bravais lattices yield no basis
if crystal_args{1} == "CUB-cP-Oh" || crystal_args{1} == "HEX-hP-D6h" || crystal_args{1} == "RHL-hR-D3d" ||...
        crystal_args{1} == "TET-tP-D4h" || crystal_args{1} == "ORC-oP-D2h" || crystal_args{1} == "MCL-mP-C2h" ||...
        crystal_args{1} == "TRI-aP-Ci"
    Gb0 = [0, 0, 0];
% - B.2 - The 2 Base-Centered Bravais Lattices turn into base-centered in reciprocal space
elseif crystal_args{1} == "ORC-oS-D2h" || crystal_args{1} == "MCL-mS-C2h"
    Gb0 = [0, 0, 0;...
        0.5*(G1+G2)];
% - B.3 - The 3 Body-Centered Bravais Lattice turns to FCC in reciprocal space
elseif crystal_args{1} == "CUB-cI-Oh" || crystal_args{1} == "TET-tI-D4h" || crystal_args{1} == "ORC-oI-D2h"
    Gb0 = [0, 0, 0;...
        0.5*(G1+G2);...
        0.5*(G1+G3);...
        0.5*(G2+G3)];
% - B.4 - The 2 Face-Centered Bravais Lattices turns to BCC in reciprocal space
elseif crystal_args{1} == "CUB-cF-Oh" || crystal_args{1} == "ORC-oF-D2h"
    Gb0 = [0, 0, 0;...
        0.5*(G1+G2+G3)];
end
%% C - Extracting the reciprocal-space lattice over multiple translations of G
l=1;
for i=-1:1
    for j=0:2
        for k=-1:1
            for n = 1:size(Gb0,1)
                G(l,:)=i*G1+j*G2+k*G3+Gb0(n,:);
                l=l+1;
            end
        end
    end
end
% - C.1 - Finding the axes limited for a 2x2 super structure as in G
if crystal == "CUB-cP-Oh" || crystal == "HEX-hP-D6h" || crystal == "RHL-hR-D3d" ||...
        crystal == "TET-tP-D4h" || crystal == "ORC-oP-D2h" || crystal == "MCL-mP-C2h" ||...
        crystal == "TRI-aP-Ci"
    Glims = 1.04*[min(G(:,1)), max(G(:,1)), min(G(:,2)), max(G(:,2)), min(G(:,3)), max(G(:,3))];
elseif crystal == "ORC-oS-D2h" || crystal == "MCL-mS-C2h" || crystal == "CUB-cI-Oh" ||...
        crystal == "TET-tI-D4h" || crystal == "ORC-oI-D2h" || crystal == "CUB-cF-Oh" ||...
        crystal == "ORC-oF-D2h"
    Glims = 1.04*[-1*norm(G1), norm(G1), 0, 2*norm(G2),-1*norm(G3), norm(G3)];
end

%% D - Finding all vector paths to link all nearest neighbour points
Gnn = [];
% - Filing through all possible points given in the G-matrix
for i = 1:length(G)-1
    waitbar(i/length(G)-1, wbar, 'Extracting reciprocal-space lattice...', 'Name', 'extract_lattice');
    for j = i+1:length(G)
        % -- Finding the Euclidean distance between the two lattice points
        dist  = sqrt((G(i,1)-G(j,1)).^2+(G(i,2)-G(j,2)).^2+(G(i,3)-G(j,3)).^2);
        %% 4.1 - Applying geometrical constraints for nearest neighbour points

        % --- For the 7 Primitive (Simple) Reciprocal lattices
        % Points are linked if their distance is identical to the
        % length of G1, G2 or G3.
        if crystal_args{1} == "CUB-cP-Oh" || crystal_args{1} == "HEX-hP-D6h" || crystal_args{1} == "RHL-hR-D3d" ||...
                crystal_args{1} == "TET-tP-D4h" || crystal_args{1} == "ORC-oP-D2h" || crystal_args{1} == "MCL-mP-C2h" ||...
                crystal_args{1} == "TRI-aP-Ci"
            if (abs(dist-norm(G1))) <= 0.0001 || (abs(dist-norm(G2))) <= 0.0001 || (abs(dist-norm(G3))) <= 0.0001                
                Gnn = vertcat(Gnn, [[G(i,1);G(j,1)], [G(i,2);G(j,2)], [G(i,3);G(j,3)]]);
            end

        % --- For the 2 Base-Centered Reciprocal Lattices
        % Points are linked if their distance is equal to the length
        % of the BCC basis vector. Also, if their length is equal to 
        % G1, G2 or G3, with the condition that it is not divisible by
        % the BCC basis vector (as the base-centered points are only
        % joined within the horizontal plane of the vertices).
        elseif crystal_args{1} == "ORC-oS-D2h" || crystal_args{1} == "MCL-mS-C2h" 
            if abs(dist-norm(Gb0(2,:))) <= 0.0001
                Gnn = vertcat(Gnn, [[G(i,1);G(j,1)], [G(i,2);G(j,2)], [G(i,3);G(j,3)]]);
            elseif (mod(G(i,2), 2*Gb0(2,2)) == 0 ||  abs(G(i,2)) <= 0.0001 ) && ((abs(dist-norm(G1))) <= 0.0001 || (abs(dist-norm(G2))) <= 0.0001 || (abs(dist-norm(G3))) <= 0.0001)
                Gnn = vertcat(Gnn, [[G(i,1);G(j,1)], [G(i,2);G(j,2)], [G(i,3);G(j,3)]]);
            end

        % --- For the 2 Body-Centered Reciprocal Lattices
        % Points are linked if their distance is equal to the length
        % of the BCC basis vector. Also, if their length is equal to 
        % G1, G2 or G3, with the condition that it is not divisible by
        % the BCC basis vector (as the body-centered points are only
        % joined radially outwards to vertices of unit cell).
        elseif crystal_args{1} == "CUB-cF-Oh" || crystal_args{1} == "ORC-oF-D2h"
            if abs(dist-norm(Gb0(2,:))) <= 0.0001
                Gnn = vertcat(Gnn, [[G(i,1);G(j,1)], [G(i,2);G(j,2)], [G(i,3);G(j,3)]]);
            elseif (mod(G(i,3), 2*Gb0(2,3)) == 0 || abs(G(i,3)) <= 0.001) && ((abs(dist-norm(G1))) <= 0.0001 || (abs(dist-norm(G2))) <= 0.0001 || (abs(dist-norm(G3))) <= 0.0001)
                Gnn = vertcat(Gnn, [[G(i,1);G(j,1)], [G(i,2);G(j,2)], [G(i,3);G(j,3)]]);
            elseif (mod(G(i,2), 2*Gb0(2,2)) == 0 || abs(G(i,2)) <= 0.001) && ((abs(dist-norm(G1))) <= 0.0001 || (abs(dist-norm(G2))) <= 0.0001 || (abs(dist-norm(G3))) <= 0.0001)
                Gnn = vertcat(Gnn, [[G(i,1);G(j,1)], [G(i,2);G(j,2)], [G(i,3);G(j,3)]]);
            end

        % --- For the 3 Face-Centered Reciprocal Lattices
        % Points are linked if their distance is equal to the length
        % of the FCC basis vector, with the proviso that the two points
        % are coplanar and not an FCC point. This is applied along each
        % dimension. Furthermore, if their length is equal to G1, G2 or G3, 
        % with the condition that it is not divisible by the FCC basis vector
        % (as the face-centered points are only joined radially outwards to
        % vertices of unit cell, not to other FCC points).
        elseif crystal_args{1} == "CUB-cI-Oh" || crystal_args{1} == "TET-tI-D4h" || crystal_args{1} == "ORC-oI-D2h"
            if abs(dist-norm(Gb0(2,:))) <= 0.001 && abs(G(i,1)-G(j,1)) <= 0.001 && (abs(G(i,1)) <= 0.001 || mod(G(i,1), 2*Gb0(2,1)) == 0)
                Gnn = vertcat(Gnn, [[G(i,1);G(j,1)], [G(i,2);G(j,2)], [G(i,3);G(j,3)]]);
            elseif abs(dist-norm(Gb0(2,:))) <= 0.001 && abs(G(i,2)-G(j,2)) <= 0.001 && (abs(G(i,2)) <= 0.001 || mod(G(i,2), 2*Gb0(2,2)) == 0)
                Gnn = vertcat(Gnn, [[G(i,1);G(j,1)], [G(i,2);G(j,2)], [G(i,3);G(j,3)]]);
            elseif abs(dist-norm(Gb0(3,:))) <= 0.001 && abs(G(i,3)-G(j,3)) <= 0.001 && (abs(G(i,3)) <= 0.001 || mod(G(i,3), 2*Gb0(3,3)) == 0)
                Gnn = vertcat(Gnn, [[G(i,1);G(j,1)], [G(i,2);G(j,2)], [G(i,3);G(j,3)]]);
            elseif (abs(dist-norm(G1))) <= 0.0001 && abs(G(i,1)-G(j,1)) <= 0.001 && (abs(G(i,1)) <= 0.001 || mod(G(i,1), 2*Gb0(2,1)) == 0) && (abs(G(i,3)) <= 0.001 || mod(G(i,3), 2*Gb0(3,3)) == 0)
                Gnn = vertcat(Gnn, [[G(i,1);G(j,1)], [G(i,2);G(j,2)], [G(i,3);G(j,3)]]);
            elseif (abs(dist-norm(G2))) <= 0.0001 && abs(G(i,2)-G(j,2)) <= 0.001 && (abs(G(i,2)) <= 0.001 || mod(G(i,2), 2*Gb0(2,2)) == 0) && (abs(G(i,3)) <= 0.001 || mod(G(i,3), 2*Gb0(3,3)) == 0)
                Gnn = vertcat(Gnn, [[G(i,1);G(j,1)], [G(i,2);G(j,2)], [G(i,3);G(j,3)]]);
            end
        end
    end
end
%% E - Computing the Voronoi diagram, which is identical to the First Brilluoin Zone
% - Calculating Voronoi diagram
[G_c, G_v] = voronoin(G); 
G_nx = G_c(G_v{floor(l/2)},:);
G_cell = convhulln(G_nx);
% - Defining the G_bz variable and shifting Brilluoin zone to +G2
for i = 1:size(G_cell,1)
    Gbz{i,1} = G_nx(G_cell(i,:),1)-G(floor(l/2),1)+G2(1);
    Gbz{i,2} = G_nx(G_cell(i,:),2)-G(floor(l/2),2)+G2(2);
    Gbz{i,3} = G_nx(G_cell(i,:),3)-G(floor(l/2),3)+G2(3);
end

%% Assigning the real- and reciprocal-information to MATLAB data structure
waitbar(1, wbar, 'Assigning real- and reciprocal-variables...', 'Name', 'extract_lattice');
% - Real-space information
realStr.crystal = crystal;
realStr.T1 = T1;
realStr.T2 = T2;
realStr.T3 = T3;
realStr.Tb0 = Tb0;
realStr.Tvol = Tvol;
realStr.T = T;
realStr.Tlims = Tlims;
realStr.Tnn = Tnn;
realStr.Tws = Tws;
% - Reciprocal-space information
reciStr.crystal = crystal;
reciStr.G1 = G1;
reciStr.G2 = G2;
reciStr.G3 = G3;
reciStr.Gb0 = Gb0;
reciStr.Gvol = Gvol;
reciStr.G = G;
reciStr.Glims = Glims;
reciStr.Gnn = Gnn;
reciStr.Gbz = Gbz;
%% Close wait-bar
close(wbar);

%% Figure summary of the real- and reciprocal-space structures
if plot_results == 1
    fig = figure(); set(fig, 'position', [1,1,950,650]);

    % - REAL-SPACE FIGURE
    subplot(1,2,1); hold on;
    % -- Plotting the axes lines
    line([min(T(:,1)), max(T(:,1))], [0 0], [0,0],  'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
    line([0 0], [min(T(:,2)), max(T(:,2))], [0, 0], 'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
    line([0 0], [0, 0], [min(T(:,3)), max(T(:,3))], 'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
    % -- Plotting the translational real-space lattice vectors
    plot3([0,T1(1)],[0,T1(2)],[0,T1(3)], '-', 'linewidth', 3, 'color', [0.8, 0.3, 0.3]);
    plot3([0,T2(1)],[0,T2(2)],[0,T2(3)], '-', 'linewidth', 3, 'color', [0.8, 0.3, 0.3]);
    plot3([0,T3(1)],[0,T3(2)],[0,T3(3)], '-', 'linewidth', 3, 'color', [0.8, 0.3, 0.3]);
    % -- Plotting the multiple translations of the real-space vectors
    plot3(T(:,1),T(:,2),T(:,3),'b.','markersize',20);
    % -- Plotting all the nearest neighbour points
    for i = 1:2:size(Tnn,1)
        plot3([Tnn(i,1), Tnn(i+1,1)], [Tnn(i,2), Tnn(i+1,2)], [Tnn(i,3), Tnn(i+1,3)], 'k-', 'linewidth', 1);
    end
    % -- Plotting the Wigner-Seitz cell
    for i = 1:size(Tws,1)
        patch(Tws{i,1}, Tws{i,2}, Tws{i,3},[0.8 0.3 0.3], 'FaceAlpha',0.75, 'EdgeColor', 'none');
    end
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
    xlabel('$$ \bf  x (\AA) $$', 'Interpreter', 'latex');
    ylabel('$$ \bf  y (\AA) $$', 'Interpreter', 'latex');
    zlabel('$$ \bf  z (\AA) $$', 'Interpreter', 'latex');
    xticks(round(-10*norm(T1):norm(T1):10*norm(T1),2));
    yticks(round(-10*norm(T2):norm(T2):10*norm(T2),2));
    zticks(round(-10*norm(T3):norm(T3):10*norm(T3),2));
    % -- Modify the view
    view(3); camlight(-30,24);
    title_txt1 = sprintf("%s; Real; Wigner-Seitz", crystal);
    title(title_txt1);
    axis tight equal vis3d; rotate3d on;
    pbaspect([1,1,1]);
    axis(Tlims);

    % - RECIPROCAL-SPACE FIGURE
    subplot(1,2,2); hold on;
    % -- Plotting the axes lines
    line([min(G(:,1)), max(G(:,1))], [0 0], [0,0],  'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
    line([0 0], [min(G(:,2)), max(G(:,2))], [0, 0], 'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
    line([0 0], [0, 0], [min(G(:,3)), max(G(:,3))], 'Color', [0 0 0], 'LineWidth', 0.75, 'Linestyle', '--');
    % -- Plotting the translational real-space lattice vectors
    plot3([0,G1(1)],[0,G1(2)],[0,G1(3)], '-', 'linewidth', 3, 'color', [0.3, 0.8, 0.3]);
    plot3([0,G2(1)],[0,G2(2)],[0,G2(3)], '-', 'linewidth', 3, 'color', [0.3, 0.8, 0.3]);
    plot3([0,G3(1)],[0,G3(2)],[0,G3(3)], '-', 'linewidth', 3, 'color', [0.3, 0.8, 0.3]);
    % -- Plotting the multiple translations of the real-space vectors
    plot3(G(:,1),G(:,2),G(:,3),'b.','markersize',20);
    % -- Plotting all the nearest neighbour points
    for i = 1:2:size(Gnn,1)
        plot3([Gnn(i,1), Gnn(i+1,1)], [Gnn(i,2), Gnn(i+1,2)], [Gnn(i,3), Gnn(i+1,3)], 'k-', 'linewidth', 1);
    end
    % -- Plotting the first Brilluoin Zone
    for i = 1:size(Gbz,1)
        patch(Gbz{i,1}, Gbz{i,2}, Gbz{i,3},[0.3 0.8 0.3], 'FaceAlpha',0.75, 'EdgeColor', 'none');
    end
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
    axis tight  equal vis3d; rotate3d on;
    pbaspect([1,1,1]);
    axis(Glims);
end
end