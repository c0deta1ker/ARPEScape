function [XMIN, XMAX] = ExLoc(X,Y)
% [XMIN, XMAX] = ExLoc(X,Y) finds all local extrema and interpolates their position
% inputs/outputs:
%   X,Y - vectors of x- and y-values
%   XMIN, XMAX - interpolated x-values of the local minima and maxima

% size check and reshape to rows 
nx=length(X); ny=length(Y); if (nx-ny)~=0 error('Mismatch of input vectors'); end
X=reshape(X,[1 nx]);  Y=reshape(Y,[1 ny]);

% - locating the extrema
L=Y(1:size(Y,2)-2); C=Y(2:size(Y,2)-1); R=Y(3:size(Y,2));
IMIN=1+find(C<L&C<R); IMAX=1+find(C>L&C>R);

% - interpolating
XMIN=[];
for imin=IMIN
   [x dummy]=XYMax( [Y(imin-1) Y(imin) Y(imin+1)] );
   XMIN=[XMIN X(imin-1)+x*(X(imin)-X(imin-1))];   
end
XMAX=[];
for imax=IMAX
   [x dummy]=XYMax( [Y(imax-1) Y(imax) Y(imax+1)] ); 
   XMAX=[XMAX X(imax-1)+x*(X(imax)-X(imax-1))];
end