function [fit1,Ef] = find_Ef_ADRESS(Angle,Energy,Data,Angle_window,E_window)
% [fit1,Ef] = find_Ef_ADRESS(Angle,Energy,Data,Angle_window, E_window, varargin)
% This function fits a Fermi-Dirac distribution function, multiplied by a 
% linear function, to the Data. Energy range [Energy_low, Energy_high].
% Angle range [Angle_low, Angle_high].
% This function performs a single fit, after binning all the channels.
% Note: the length of E or graph (Energy scale) can take arbitrary values
%       Channel range can also be set to 1. 

%% Default parameters
% -- Estimating the energy window if no entry is entered
if nargin < 5; E_window = []; end
if isempty(E_window)
    data                = squeeze(Data);
    [~, EF, ~]          = AlignEF(data, Energy, 0);
    E_window            = [EF-0.6 EF+0.6];
end 
% -- Estimating the angle window if no entry is entered
if nargin < 4; Angle_window = []; end
if isempty(Angle_window)
    Angle_window = [min(Angle(:))+0.25, max(Angle(:))-0.25];
end 
%% - 1 - Initialising variables
% Extracting the angle window
Angle_low  	= min(Angle_window);
Angle_high  = max(Angle_window);
if size(Angle, 1) ~= 1; Data = Data(Angle>Angle_low&Angle<Angle_high); 
else; Data = Data(:,Angle>Angle_low&Angle<Angle_high); 
end

%% - 2 - Extracting the index of energy window
[~,index_min]   = min(abs(Energy-min(E_window)));
[~,index_max]   = min(abs(Energy-max(E_window)));
channel_start   = index_min;
channel_end     = index_max;
E               = Energy;
graph           = Data;
%% - 3 - Summing the data over all the bins of the energy window
graph_sum   = sum(graph(channel_start:channel_end,:),2);
graph_sum   = graph_sum / max(graph_sum);
E           = E(channel_start:channel_end);
% -- Ensuring both E and graph_sum are of the same size
if size(E, 2) >1; E = E'; end
if size(graph_sum, 2) >1; graph_sum = graph_sum'; end

%% - 4 - Executing the Fermi fit
% Defining the function
fermi_fit       = fittype(@(e, a, b, c, d, x) e+(a*(x-b)+d)./(1+exp((x-b)/c)));
% Initial fit parameters
dgraph_sum      = graph_sum(10:end)-graph_sum(1:end-9);
dgraph_sum      = [graph_sum(1)*ones(1,9)   dgraph_sum'];
Ef_guess_loc   	= find(dgraph_sum == min(dgraph_sum));
Ef_guess_loc    = Ef_guess_loc(1);
a_start         = (graph_sum(Ef_guess_loc) - graph_sum(1))/(E(Ef_guess_loc) - E(1));
b_start         = E(Ef_guess_loc);
c_start         = .1;
d_start         = max(graph_sum);
e_start         = min(graph_sum);
% Executing the fitting operation
[fit1,gof1,out1]=fit(E,graph_sum,fermi_fit,'start',[e_start a_start b_start c_start d_start]);
fit1para = coeffvalues(fit1);
Ef = fit1para(3);
yy=feval(fermi_fit,fit1para(1),fit1para(2),fit1para(3),fit1para(4),fit1para(5),E);
%% - 5 - Plotting the best fit FDD
fig = figure(); 
fig.Position(3) = 1.0*400;
fig.Position(4) = 1.0*400;
hold on;
plot(E, graph_sum,'o', 'Color', 'r', 'markerfacecolor', 'r');
plot(E, yy, 'k-', 'linewidth', 2.5, 'linestyle', '-');
line([1, 1].*Ef, [min(yy(:)), max(yy(:))], 'color', 'b', 'linewidth', 1.5, 'linestyle', ':'); 
% -- Formatting the figure
% axis padded; 
gca_props(); grid on;
ylabel('$$ \bf Intensity $$', 'Interpreter', 'latex');
xlabel('$$ \bf E_B [eV] $$', 'Interpreter', 'latex');
legend({'Data', 'Fit', 'Ef'}, 'location', 'best', 'interpreter', 'none');
