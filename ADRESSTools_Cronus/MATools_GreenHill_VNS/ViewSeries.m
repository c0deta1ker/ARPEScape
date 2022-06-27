function ViewSeries(ACorr,EAlign,Scan,Data,nRows,nCols,colorScheme,minFrac,maxFrac)
% View series of downsampled images. The inputs are ACorr (2D or 3D array of [warping 
% corrected] angles), EAlign (2D or 3D matrix of curvature corrected/aligned to EF energies), 
% Scan (row), Data (3D data array Data), nRows and nCols as number of rows and columns 
% of the frames in which the images are shown, and their colorScheme from the standard 
% Matlab colormaps including 'jet', 'hot' and 'gray'. minFrac and maxFrac are contrast enhancement 
% parameters used by SetContrast.m inside.  

% downsampling parameters
[nR,nC,~]=size(EAlign); nR=nRows*floor(nR/nRows); nC=nCols*floor(nC/nCols);

% contrast parameters
if nargin<8; minFrac=1e-2; maxFrac=1-1e-4; end

% frame loop
for ipage=1:ceil(size(Data,3)/nRows/nCols)
   figure; orient tall 
   for iframe=1:nRows*nCols
      n=(ipage-1)*nRows*nCols+iframe; if n>size(Scan,2); break; end
% frames
      if ndims(ACorr)<3; ACorrFrame=ACorr; else ACorrFrame=squeeze(ACorr(:,:,n)); end
      if ndims(EAlign)<3; EAlignFrame=EAlign; else EAlignFrame=squeeze(EAlign(:,:,n)); end
      DataFrame=squeeze(Data(:,:,n));
% downsampling
% - angles
      ACorrFrame1=0; for iF=1:nCols; ACorrFrame1=ACorrFrame1+ACorrFrame(:,iF:nCols:nC)/nCols; end
      ACorrFrame=0; for iF=1:nRows; ACorrFrame=ACorrFrame+ACorrFrame1(iF:nRows:nR,:)/nRows; end
% - energies      
      EAlignFrame1=0; for iF=1:nCols; EAlignFrame1=EAlignFrame1+EAlignFrame(:,iF:nCols:nC)/nCols; end
      EAlignFrame=0; for iF=1:nRows; EAlignFrame=EAlignFrame+EAlignFrame1(iF:nRows:nR,:)/nRows; end
% - data      
      DataFrame1=0; for iF=2:nCols; DataFrame1=DataFrame1+DataFrame(:,iF:nCols:nC)/nCols; end
      DataFrame=0; for iF=2:nRows; DataFrame=DataFrame+DataFrame1(iF:nRows:nR,:)/nRows; end
% contrast enhancement 
      DataFrame=SetContrast(DataFrame,minFrac,maxFrac);
% plot
      subplot(nRows,nCols,iframe); 
%      h=pcolor(ACorrFrame,EAlignFrame,DataFrame); set(h,'EdgeColor','None','FaceColor','Interp'); 
      ImData(ACorrFrame,EAlignFrame,DataFrame);
      colormap('jet'); if nargin>6; eval(['colormap ''' colorScheme '''']); end
      title(['Scan=' num2str(Scan(n))],'Interpreter','none','fontweight','normal'); 
   end
end