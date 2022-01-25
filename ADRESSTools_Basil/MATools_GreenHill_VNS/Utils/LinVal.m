function Z=LinVal(FString,Coeffs,X,Y)
% Z=LinVal(FString,Coeffs,X,Y) returns values of the linear combination of 
% the functions of X and Y defined in FString (space separated) with 
% the coefficients defined in Coeffs. 
% Example: Z=LinVal('1 X.^2 X.Y Y.^2',Coeffs,X,Y)

% squeeze multiple spaces from FString
while ~isempty(strfind(FString,'  '))
   FString=strrep(FString,'  ',' ');
end

% consistency check
if length(strfind(FString,' '))~=length(Coeffs)-1
    Z=[]; disp('LinVal error: FString misfits with Coeffs'); return
end

% formation of the linear combination string
FString=['Coeffs(1)*' FString]; % leading term
for iCoeff=2:length(Coeffs)
   iSep=strfind(FString,' '); iSep=iSep(1);
   FString=[FString(1:iSep-1) '+Coeffs(' num2str(iCoeff) ')*' FString(iSep+1:end)];
end
    
% values
eval(['Z=' FString ';']);