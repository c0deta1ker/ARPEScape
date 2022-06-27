function D=DifC(X,n)
% D = DifC(X,n) calculates the central differencies for the vector X of 
% the n=1 or 2 order. The result is end-padded to keep the length of X. 
% If A is matrix, the differences are calculated along the columns.
 
% consistency check
if n~=1 && n~=2; disp('DifC: Error, only n=1 or 2 are supported'); return; end    

% transpose the row vector input
if size(X,1)==1; X=X'; tra=1; else tra=0; end

% central differences
if n==1
   D=( X(3:end,:)-X(1:end-2,:) )/2;
else
   D=X(3:end,:)+X(1:end-2,:)-2*X(2:end-1,:);
end
% padding
D=[D(1,:); D; D(end,:)];

% back transpose when the row vector input
if tra==1; D=D'; end