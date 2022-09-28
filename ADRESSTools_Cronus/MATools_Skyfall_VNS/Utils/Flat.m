function [AFlat,Coeffs]=Flat(Angle,Scan,Slice)
% AFlat=Flat(Angle,Scan,Slice) returns a bilinear fit of Slice

% flatten the arrays
SliceF=Slice(:);
[AngleM,ScanM]=meshgrid(Angle,Scan); AngleF=AngleM(:); ScanF=ScanM(:); 
OnesF=ones(size(SliceF));
% fitting
Coeffs=[AngleF ScanF OnesF]\SliceF;
% fitting plane
AFlat=Coeffs(1)*AngleM+Coeffs(2)*ScanM+Coeffs(3);
% figure; imagesc(Angle,Scan,Slice); set(gca,'TickDir','Out'); axis equal tight;
% figure; imagesc(Angle,Scan,FitM); set(gca,'TickDir','Out'); axis equal tight;

