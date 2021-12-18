function Fit=LinFit(X,Y,Z,FString)
% linear fit of Z(X,Y) with the functions of X and Y specified in FString 
% (space separated). The constant term is represented by leading 1.
% Example: Fit=LinFit(X,Y,Z,'1 X.^2 X.*Y Y.^2')

% replacing the constant term
if isequal(FString(1),'1'); FString=['ones(size(X))' FString(2:end)]; end

% flattening
X=X(:); Y=Y(:); Z=Z(:); 
if length(Y)~=length(X) || length(Z)~=length(X) 
   Fit=[]; disp('LinFit error: Incompatible input arrays'); return
end

% - fitting
eval(['Fit=[' FString ']\Z;']); Fit=Fit';