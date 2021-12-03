function [Angle_axis_truncated, Energy_axis_truncated, Data_truncated]=...
    truncate_ADRESS_data(Angle, Energy, Data, Angle_window, Energy_window)
% truncates data and axis in range given by windows
% [Angle_axis_truncated, Energy_axis_truncated, Data_truncated]=truncate_ADRESS_data(Angle, Energy, Data, Angle_window, Energy_window)
Angle_low=min(Angle_window);
Angle_high=max(Angle_window);
E_low=min(Energy_window);
E_high=max(Energy_window);

Energy_axis_truncated=Energy(Energy>E_low&Energy<E_high); 
Angle_axis_truncated=Angle(Angle>Angle_low&Angle<Angle_high);
if ndims(Data)>2
    Data_truncated=Data(Energy>E_low&Energy<E_high,Angle>Angle_low&Angle<Angle_high,:); 
else
    Data_truncated=Data(Energy>E_low&Energy<E_high,Angle>Angle_low&Angle<Angle_high); 
end    