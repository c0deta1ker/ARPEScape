function h=ImData(X,Y,Z,style)
% [h=] ImData(X,Y,Z[,style]) plots image of Z(X,Y) using the style 'flat' (default) 
% or 'interp' style. X and Y are vectors or 2D arrays, Z is 2D array.
% h (optional) can be used for further handling of the image (rotation, etc) 

% default style
if nargin==3; style='flat'; end

% squeeze input arrays to 2D
% X=squeeze(X); Y=squeeze(X); Z=squeeze(Z); 

% X and Y matrices
% [X,Y]=meshgrid(X,Y); end
if size(X,1)==1||size(X,2)==1 
   if size(X,1)>size(X,2); X=X'; end  
   X=repmat(X,size(Z,1),1);
end
if size(Y,1)==1||size(Y,2)==1
   if size(Y,1)<size(Y,2); Y=Y'; end 
   Y=repmat(Y,1,size(Z,2));
end

if sum(size(Z)-size(X))~=0 || sum(size(Z)-size(Y))~=0 error('Mismatch of input matrices'); end 
    
switch lower(style)
case('flat')    
% appending along X
   DX=diff(X,1,2)/2; XV=[X(:,1)-DX(:,1) X(:,1:end-1)+DX X(:,end)+DX(:,end)];
   DXV=XV(end,:)-XV(end-1,:); XV=[XV; XV(end,:)+DXV];
% appending along Y
   DY=diff(Y,1,1)/2; YV=[Y(1,:)-DY(1,:); Y(1:end-1,:)+DY; Y(end,:)+DY(end,:)];
   DYV=YV(:,end)-YV(:,end-1); YV=[YV YV(:,end)+DYV];
% appending Z
   ZV=[Z NaN*ones(size(Z,1),1)]; ZV=[ZV; NaN*ones(1,size(ZV,2))];
% plotting
   h=pcolor(XV,YV,ZV); set(h,'EdgeColor','None','FaceColor','Flat');
case ('interp')
   h=pcolor(X,Y,Z); set(h,'EdgeColor','None','FaceColor','Interp');
otherwise
   error('Invalid image style')
end

% axes on top and ticks
set(gca,'Layer','Top','TickDir','Out');