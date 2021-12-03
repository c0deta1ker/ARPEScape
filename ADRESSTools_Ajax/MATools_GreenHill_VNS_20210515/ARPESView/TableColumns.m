% X. Wang, 02-Feb-2021
classdef TableColumns < uint32
    enumeration 
        % File
        File(1)
        Comment(2)
        % Beamline
        hv(3)
        Polarization(4)
        Slit(5)
        % Analyzer
        Mode(6)
        Epass(7)
        EbEk(8)
        Energy(9)
        dt(10)
        Sweeps(11)
        % Manipulator
        Theta(12)
        Tilt(13)
        Azimuth(14)
        X(15)
        Y(16)
        Z(17)
        % Control
        Activate(18)
        Action(19)
        Time(20)
    end
    methods(Static)
        function v = numberColumns()
            v = 20;
        end
        function v = vectorizableColumns()
            v = [hv, Theat, Tilt, Azimuth, X, Z];
        end
        function v = emptyRow()
            v = {'' '' ...              % File
                '' '' [] ...            % Beamline
                '' [] '' '' '' [] ...   % Analyzer
                '' '' '' '' '' ''...    % Manipulator
                false 'Select' ''      % Control
            };
        end
        function v = columnNames()
            v = {'File', 'Comment', ...
                '<html><b><i>hv</i></b></html>', 'Pol', 'Slit', ...
                'Mode', 'Epass', 'Eb/k', 'Energy', 'dt(s)', 'Sweeps', ...
                '<html><b><i>Theta</i></b></html>', ...
                '<html><b><i>Tilt</i></b></html>', ...
                '<html><b><i>Azim</i></b></html>', ...
                '<html><b><i>X</i></b></html>', ...
                '<html><b><i>Y</i></b></html>', ...
                '<html><b><i>Z</i></b></html>', ...
                'Activate', 'Action', 'T(min)'};
        end
        function v = columnFormats()
            v = {'char', 'char', ...
                'char', {'LH','LV','C+','C-'}, 'char',...
                {'WAM','MAM','LAD'}, 'char', {'Eb','Ek'}, 'char', 'char', 'char', ...
                'char', 'char', 'char', 'char', 'char', 'char' ...
                'logical', {'Select','Clear','Copy','Paste','Edit Raw', 'Load Raw'}, 'char'};
        end
        function v = columnWidth()
            v = {122 122 ...
                80 40 40 ...
                45 40 40 80 30 0 ...
                70 70 70 70 0 70 ...
                50 70 70};
        end
    end
end