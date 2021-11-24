function [XMin,YMin,XMax,YMax] = ExtremeXY(X,Y)
% [XMin,YMin,XMax,YMax] = ExtremeXY(X,Y) finds and interpolates all local extremes
% Inputs/outputs:
%   X,Y - vectors (rows or columns) of x- and y-values
%   XMin,YMin,XMax,YMax - interpolated x- and y-values of the local minima and maxima
% Note: X shall not necessarily be equidistant, but plaid and not repeating values 

% locating indices of the extremes
L=Y(1:end-2); C=Y(2:end-1); R=Y(3:end);
IMin=1+find(C<L&C<R); IMax=1+find(C>L&C>R);

% offsets to pre-condition the polynom fit
xOffset=mean(X); yOffset=mean(Y);
X=X-xOffset; Y=Y-yOffset;

% interpolating minima and offset backcorrection
XMin=[]; YMin=[];
if ~isempty(IMin) 
   for iMin=IMin
      X(iMin-1:iMin+1); Y(iMin-1:iMin+1);
      Fit=polyfit(X(iMin-1:iMin+1),Y(iMin-1:iMin+1),2);
      xMin=-0.5*Fit(2)/Fit(1); XMin=[XMin xMin];
      yMin=polyval(Fit,xMin); YMin=[YMin yMin];
   end; 
   XMin=XMin+xOffset; YMin=YMin+yOffset;
end

% interpolating maxima and offset backcorrection
XMax=[]; YMax=[];
if ~isempty(IMax)
   for iMax=IMax
      Fit=polyfit(X(iMax-1:iMax+1),Y(iMax-1:iMax+1),2);
      xMax=-0.5*Fit(2)/Fit(1); XMax=[XMax xMax];
      yMax=polyval(Fit,xMax); YMax=[YMax yMax];
   end; 
   XMax=XMax+xOffset; YMax=YMax+yOffset;
end
 