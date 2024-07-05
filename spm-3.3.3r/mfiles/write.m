function write(SPM)

%
% input write function for iSOSIA 3
%
% DLE 19/1/2015
%


%File mesh.input
fid = fopen('./input/mesh.input','w');
fprintf(fid,'%5.7e\n',SPM.mesh.L);
fprintf(fid,'%5.7e\n',SPM.mesh.H);
fprintf(fid,'%d\n',SPM.mesh.nx);
fprintf(fid,'%d\n',SPM.mesh.ny);
fprintf(fid,'%5.7e\n',SPM.mesh.dx);
fprintf(fid,'%5.7e\n',SPM.mesh.dy);
fprintf(fid,'%5.7e\n',SPM.mesh.hmin);
fprintf(fid,'%5.7e\n',SPM.mesh.ct);
fprintf(fid,'%5.7e\n',SPM.mesh.maxtime);
fprintf(fid,'%5.7e\n',SPM.mesh.filetime);
fprintf(fid,'%5.7e\n',SPM.mesh.maxdt);
fprintf(fid,'%5.7e\n',SPM.mesh.gravity);
fprintf(fid,'%5.7e\n',SPM.mesh.maxb);
fprintf(fid,'%5.7e\n',SPM.mesh.maxs);
fprintf(fid,'%d\n',SPM.mesh.nmoni);
fprintf(fid,'%d\n',SPM.mesh.periodic);
fprintf(fid,'%d\n',SPM.mesh.slidingmode);
fprintf(fid,'%d\n',SPM.mesh.coldbased);
fprintf(fid,'%d\n',SPM.mesh.hydromode);
fprintf(fid,'%d\n',SPM.mesh.dofluvial);
fprintf(fid,'%d\n',SPM.mesh.doice);
fprintf(fid,'%d\n',SPM.mesh.dosliding);
fprintf(fid,'%d\n',SPM.mesh.doperiglacial);
fprintf(fid,'%d\n',SPM.mesh.doglacialhydrology);
fprintf(fid,'%d\n',SPM.mesh.doglacialerosion);
fprintf(fid,'%d\n',SPM.mesh.dohillslope);
fprintf(fid,'%d\n',SPM.mesh.dohillslopeerosion);
fprintf(fid,'%d\n',SPM.mesh.dolandslide);
fprintf(fid,'%d\n',SPM.mesh.doavalance);
fprintf(fid,'%d\n',SPM.mesh.doisostasy);
fprintf(fid,'%d\n',SPM.mesh.docelldata);
fprintf(fid,'%d\n',SPM.mesh.dosediment);
fprintf(fid,'%d\n',SPM.mesh.doglacialsedi);
fprintf(fid,'%d\n',SPM.mesh.dodebrisablation);
fclose(fid);

%File iprop.input
fid = fopen('./input/iprop.input','w');
fprintf(fid,'%5.7e\n',SPM.iprop.gamma);
fprintf(fid,'%5.7e\n',SPM.iprop.gamma0);
fprintf(fid,'%5.7e\n',SPM.iprop.Cs);
fprintf(fid,'%5.7e\n',SPM.iprop.latentheat);
fprintf(fid,'%5.7e\n',SPM.iprop.ki);
fprintf(fid,'%5.7e\n',SPM.iprop.rho);
fprintf(fid,'%5.7e\n',SPM.iprop.cp);
fprintf(fid,'%5.7e\n',SPM.iprop.ifac);
fprintf(fid,'%5.7e\n',SPM.iprop.sfac);
fprintf(fid,'%5.7e\n',SPM.iprop.vbfac);
fprintf(fid,'%d\n',SPM.iprop.maxitt_v);
fprintf(fid,'%d\n',SPM.iprop.maxitt_s);
fprintf(fid,'%5.7e\n',SPM.iprop.C);
fprintf(fid,'%5.7e\n',SPM.iprop.L0);
fprintf(fid,'%5.7e\n',SPM.iprop.minbeta);
fprintf(fid,'%5.7e\n',SPM.iprop.maxsliding);
fprintf(fid,'%5.7e\n',SPM.iprop.maxdeformation);
fprintf(fid,'%5.7e\n',SPM.iprop.mf);
fprintf(fid,'%5.7e\n',SPM.iprop.qfac);
fprintf(fid,'%5.7e\n',SPM.iprop.Ka);
fprintf(fid,'%5.7e\n',SPM.iprop.ap);
fclose(fid);

%File mprop.input
fid = fopen('./input/mprop.input','wb');
fwrite(fid,SPM.mprop.mtype,'int');
fwrite(fid,SPM.mprop.avaslope,'double');
fwrite(fid,SPM.mprop.avacurv,'double');
fwrite(fid,SPM.mprop.maxslope,'double');
fwrite(fid,SPM.mprop.lrate,'double');
fwrite(fid,SPM.mprop.Tsl,'double');
fwrite(fid,SPM.mprop.maxacc,'double');
fwrite(fid,SPM.mprop.maxabla,'double');
fwrite(fid,SPM.mprop.accgrad,'double');
fwrite(fid,SPM.mprop.ablgrad,'double');
fwrite(fid,SPM.mprop.qb,'double');
fwrite(fid,SPM.mprop.sedifrac,'double');
fwrite(fid,SPM.mprop.Ldebris,'double');
fwrite(fid,length(SPM.mprop.Temp(1,:)),'int');
fwrite(fid,SPM.mprop.Temp,'double');
fwrite(fid,length(SPM.mprop.Mrate(1,:)),'int');
fwrite(fid,SPM.mprop.Mrate,'double');
fclose(fid);

%File hwprop.input
fid = fopen('./input/hwprop.input','wb');
fwrite(fid,SPM.hwprop.Kgw,'double');
fwrite(fid,SPM.hwprop.a2w,'double');
fwrite(fid,SPM.hwprop.po,'double');
fwrite(fid,SPM.hwprop.tscale,'double');
fclose(fid);

%File hprop.input
fid = fopen('./input/hprop.input','wb');
fwrite(fid,SPM.hprop.Ks,'double');
fwrite(fid,SPM.hprop.sc,'double');
fwrite(fid,SPM.hprop.Ke,'double');
fwrite(fid,SPM.hprop.Ls,'double');
fwrite(fid,SPM.hprop.gamma,'double');
fwrite(fid,SPM.hprop.Nc,'int');
fclose(fid);

%File fprop.input
fid = fopen('./input/fprop.input','wb');
fwrite(fid,SPM.fprop.pr,'double');
fwrite(fid,SPM.fprop.rho_s,'double');
fwrite(fid,SPM.fprop.Dg,'double');
fwrite(fid,SPM.fprop.kw,'double');
fwrite(fid,SPM.fprop.tau_c,'double');
fwrite(fid,SPM.fprop.Kt,'double');
fwrite(fid,SPM.fprop.Ke,'double');
fclose(fid);

%File meshdata.input
fid = fopen('./input/meshdata.input','wb');
fwrite(fid,SPM.data.bed,'double');
fwrite(fid,SPM.data.ice,'double');
fwrite(fid,SPM.data.sedi,'double');
fwrite(fid,SPM.data.maxacc,'double');
fwrite(fid,SPM.data.mrate,'double');
fwrite(fid,SPM.data.srate,'double');
fwrite(fid,SPM.data.fixflag,'int');
fwrite(fid,SPM.data.phi,'double');
fclose(fid);

%File periglacial.input
if SPM.mesh.doperiglacial,
  fid = fopen('./input/periglacial.input','wb');
  fwrite(fid,SPM.pgprop.nHs,'int');
  fwrite(fid,SPM.pgprop.nT0,'int');
  fwrite(fid,SPM.pgprop.Hsv,'double');
  fwrite(fid,SPM.pgprop.T0v,'double');
  fwrite(fid,SPM.pgprop.Ci,'double');
  fwrite(fid,SPM.pgprop.Tr,'double');
  fwrite(fid,SPM.pgprop.rho_b,'double');
  fwrite(fid,SPM.pgprop.rho_s,'double');
  fwrite(fid,SPM.pgprop.Ke,'double');
  fwrite(fid,SPM.pgprop.Kt,'double');
  fwrite(fid,SPM.pgprop.maxsedi,'double');
  fwrite(fid,SPM.pgprop.maxice,'double');
  fwrite(fid,SPM.pgprop.minslope,'double');
  fwrite(fid,SPM.pgprop.minsedi,'double');
  fclose(fid);
end;

if SPM.mesh.docelldata,
  load ./input/nodelist.mat; ncells = length(nnx);
  fid = fopen('./input/cellnumbers.input','wb');
  fwrite(fid,ncells,'int');
  for i=1:ncells,
    fwrite(fid,nnx(i),'int');
    fwrite(fid,nny(i),'int');
  end;
  fclose(fid);
end;

if SPM.mesh.doglacialsedi,
  fid = fopen('./input/Vsdata.input','wb');
  for i=1:20,
    fwrite(fid,SPM.data.Vs{i},'double');
  end;
end;

  
% write matlab structure 
save('./input/SPM.mat');

