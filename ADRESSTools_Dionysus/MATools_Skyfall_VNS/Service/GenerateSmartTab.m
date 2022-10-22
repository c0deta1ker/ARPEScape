mode='WAM'; Ek='350'; hv='1000'; pol='LV'; slitMult=1500; dT='30';

prefix=[mode 'W_Ek' Ek];
fid=fopen(['Tab' prefix '.txt'],'wt');
head='/*File|Comment|hv|Pol|Slit|Mode|Epass|Eb/Ek|Range|dt(s)|Theta|Tilt|Azimuth|X|Y|Z*/'; fprintf(fid,'%s\n',head);

%for Ep=[30:5:100 110:10:200 210:20:350] %MAM
for Ep=[30:5:100  110:10:230] % WAM
%for Ep=[30:5:100  110:10:200 210:20:350] % LAD
   slit=round(slitMult/Ep);
   line=strcat(prefix,'_Ep',num2str(Ep),'||',hv,'|',pol,'|',num2str(slit),'|',mode,'|',num2str(Ep),'|','Ek','|',Ek,'|',dT,'||||||');
   fprintf(fid,'%s\n',line);
end

fclose(fid);