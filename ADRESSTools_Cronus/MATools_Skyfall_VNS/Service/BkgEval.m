% bi-quadratic fitting of the CCD background for data acquisition

% parameters
% nCorn=32; startFit=125; endFit=10;
nCorn=32; startFit=150; endFit=5;

% read in bkg file (ReadARPES.m fails with the sp2-file not containing corrected data)
BkgFile='X:\CommiData\2014\03_2014\Bkg\SensiCamQE_Bkg.sp2';
BkgM=ReaderSP2(BkgFile); BkgM=BkgM.image;
[nr,nc]=size(BkgM);

% bkg normalization
BkgNorm=BkgM(1:nCorn,1:nCorn)+BkgM(end-nCorn+1:end,end-nCorn+1:end)+...
        BkgM(1:nCorn,end-nCorn+1:end)+BkgM(end-nCorn+1:end,1:nCorn);
BkgN=BkgM/sum(sum(BkgNorm));

% bi-quadratic fit
% - normalized arrays
EL=(0:nc-1)/(nc-1); EM=repmat(EL,nr,1);
AL=(0:nr-1)'/(nr-1); AM=repmat(AL,1,nc);
% - fitting region
ER=EM(:,startFit:end-endFit); AR=AM(:,startFit:end-endFit); BkgR=BkgN(:,startFit:end-endFit);
% - fitting
% FString='1 X X.^2 X.^3 X.*Y X.^2.*Y X.^3.*Y Y';
FString1='1 X X.^2 X.^3 X.^4 X.^5 X.^6 X.^7 X.^8 X.^9 X.^10';
FString2=' Y X.*Y X.^2.*Y X.^3.*Y X.^4.*Y X.^5.*Y X.^6.*Y X.^7.*Y X.^8.*Y X.^9.*Y';
FString3=' Y.^2 X.*Y.^2 X.^2.*Y.^2 X.^3.*Y.^2 X.^4.*Y.^2 X.^5.*Y.^2 X.^6.*Y.^2 X.^7.*Y.^2 X.^8.*Y.^2';
FString=[FString1 FString2 FString3];
Fit=LinFit(ER,AR,BkgR,FString);

% check the fit
% - input bkg
%subplot(3,1,1); plot3(AM,EM,medfilt2(BkgN,[10,10],'symmetric'))
% - fitted surface
BkgFit=LinVal(FString,Fit,EM,AM);
% - bkg fit
%subplot(3,1,2); plot3(AM,EM,BkgFit)   
% - angle-integrated bkg vs fitted surface
subplot(3,1,3); plot(1:nc,sum(BkgN,1),1:nc,sum(BkgFit))

% save the fit
save BkgFit.dat Fit -ascii

% save the model bkg in the CCDAcqure format
BkgFit=65e3*BkgFit/max(max(BkgFit));
fid=fopen('Bkg.pnm','w');
   fprintf(fid,'%s\n','P2'); 
   fprintf(fid,'%4u %4u\n%4u\n',[size(BkgFit,2) size(BkgFit,1) max(max(BkgFit))]);
   for n=1:size(BkgFit,1);
      fstring=[repmat('%4u\t ',1,size(BkgFit,2)-1) '%4u\n' ]; 
      fprintf(fid,fstring,uint16(BkgFit(n,:))); 
   end
fclose(fid);