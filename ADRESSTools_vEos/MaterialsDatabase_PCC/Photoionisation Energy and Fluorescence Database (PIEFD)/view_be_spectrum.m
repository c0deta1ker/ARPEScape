function view_be_spectrum(element)
% view_be_spectrum(element)
%   This is a function that plots the photoionisation cross-sections and
%   asymmetry parameters of all the core-levels for a particular element
%   defined by the user. This can be used to quickly view all the available 
%   data for a particular element.
%
%   IN:
%   -   element:	vector or single string of the elements; e.g. "H", "He", "Si", "In"...
%
%   OUT: (figure output)

%% Default parameters (Parameters for Silicon)
if nargin < 1; element = ["Si", "O", "Al"];  end
if isempty(element); element = ["Si", "O", "Al"]; end
% -- Ensuring the input is a string
element = convertCharsToStrings(element);

%% - Filing through all the elements and extracting binding energies
for i = 1:length(element)
    % - Extracting the PIXSA properties from the Database
    piefd_props = get_piefd_props(element(i));
    piefd_BE = table2cell(piefd_props);
    piefd_BE(1:5) = []; piefd_BE(25:end) = [];
    piefd_BE = cell2mat(piefd_BE);
    % -- Defining the x- and y-axis values
    x(:,i) = -1* piefd_BE;
    y(:,i) = 0.05*i + ones(size(x(:,i)));
end
%% - Defining the labels of the core-levels extracted
% LABEL = {...
%     'K1',...                                            % 1s
%     'L1', 'L2', 'L3',...                                % 2s, 2p
%     'M1', 'M2', 'M3', 'M4', 'M5',...                    % 3s, 3p, 3d
%     'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7',...        % 4s, 4p, 4d, 4f
%     'O1', 'O2', 'O3', 'O4', 'O5',...                    % 5s, 5p, 5d
%     'P1', 'P2', 'P3',...                                % 6s, 6p
%     };
LABEL = {...
    '1s',...                                                                        % K1
    '2s', '2p(1/2)', '2p(3/2)',...                                                  % L1,L2,L3
    '3s', '3p(1/2)', '3p(3/2)', '3d(3/2)','3d(5/2)',...                             % M1,M2,M3,M4,M5
    '4s', '4p(1/2)', '4p(3/2)', '4d(3/2)','4d(5/2)', '4f(5/2)', '4f(7/2)',...       % N1,N2,N3,N4,N5,N6,N7
    '5s', '5p(1/2)', '5p(3/2)', '5d(3/2)','5d(5/2)',...                             % O1,O2,O3,O4,O5
    '6s', '6p(1/2)', '6p(3/2)',...                                                  % P1,P2,P3
    };
yoff = linspace(0,1,length(LABEL));
%% - Plotting the binding energy spectrum
fig = figure(); hold on;
fig.Position(3) = 800; 
fig.Position(4) = 400;
% -- Plotting the binding energy spectrum
for i = 1:length(element)
    stem(x(:,i), y(:,i), '-', 'linewidth', 3.0, 'marker', 'none');
end
% -- Adding labels to each component
for i = 1:length(element)
    for j = 1:length(x(:,i))
        label = sprintf('%s(%s)', string(element(i)), string(LABEL{j}));
        text(x(j,i), y(j,i) - yoff(j), label, 'Rotation',45, 'FontWeight','bold', 'FontSize',10);
    end
end
% -- Formatting the figure
title_str = "[";
for i = 1:length(element)
    if i == length(element);    title_str = title_str + string(element(i)) + "] : BE spectrum";
    else;                       title_str = title_str + string(element(i)) + ", ";
    end
end
title(title_str);
xlabel('$$ \bf Binding\ Energy\ [eV] $$', 'interpreter', 'latex');
ylabel('$$ \bf Intensity\ [arb.] $$', 'interpreter', 'latex');
box off;
ylim([0, 1.25*max(y(:))]);

end