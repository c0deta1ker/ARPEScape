function [xMax,yMax,zMax]=MaxLoc2(X,Y,Z)
% [xMax,yMax,zMax]=MaxLoc2(X,Y,Z) finds the global maximum of Z(X,Y) and
% quadratically interpolates its coordinates and value

% locate the global maximum
ZMax=max(max(Z)); [iMax,jMax]=find(Z==ZMax);
X=X(jMax-1:jMax+1); Y=Y(iMax-1:iMax+1); Z=Z(iMax-1:iMax+1,jMax-1:jMax+1);
% test
% figure; ImData(X,Y,Z);

% quadratic fit of the global maximum 
% - flatten the arrays
[XM,YM]=meshgrid(X,Y); XF=XM(:); YF=YM(:); 
ZF=Z(:);
OnesF=ones(size(ZF));
% - fit coeffs
A=[XF.^2 YF.^2 XF.*YF XF YF OnesF]\ZF;
% - fit surface
ZFit=A(1)*XM.^2+A(2)*YM.^2+A(3)*XM.*YM+A(4)*XM+A(5)*YM+A(6);
% - output test
% figure; ImData(X,Y,ZFit);

% - coordinates and value of the global maximum
xMax=(2*A(2).*A(4)-A(3).*A(5))/(A(3).^2-4*A(1).*A(2));
yMax=(-2*A(1)*xMax-A(4))/A(3);
zMax=A(1)*xMax.^2+A(2)*yMax.^2+A(3)*xMax.*yMax+A(4)*xMax+A(5)*yMax+A(6);