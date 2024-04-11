
meshtype read_mesh(char *path,char file[200],FILE *fp)
{

  meshtype mesh;

  /*read from file mesh.input*/
  sprintf(file,"%s/input/mesh.input",path); /*filename*/
  if ((fp = fopen(file,"rt")) == NULL) {
    printf("Cannot open file : %s/input/mesh.input\n",path);
    exit(1);
  }
  fscanf(fp,"%lf\n",&mesh.L); 
  fscanf(fp,"%lf\n",&mesh.H); 
  fscanf(fp,"%d\n",&mesh.nx); 
  fscanf(fp,"%d\n",&mesh.ny); 
  fscanf(fp,"%lf\n",&mesh.dx); 
  fscanf(fp,"%lf\n",&mesh.dy); 
  fscanf(fp,"%lf\n",&mesh.hmin); 
  fscanf(fp,"%lf\n",&mesh.ct); 
  fscanf(fp,"%lf\n",&mesh.maxtime); 
  fscanf(fp,"%lf\n",&mesh.filetime); 
  fscanf(fp,"%lf\n",&mesh.maxdt); 
  fscanf(fp,"%lf\n",&mesh.gravity); 
  fscanf(fp,"%lf\n",&mesh.maxb); 
  fscanf(fp,"%lf\n",&mesh.maxs); 
  fscanf(fp,"%d\n",&mesh.nmoni); 
  fscanf(fp,"%d\n",&mesh.gmode); 
  fscanf(fp,"%d\n",&mesh.slidingmode); 
  fscanf(fp,"%d\n",&mesh.coldbased); 
  fscanf(fp,"%d\n",&mesh.hydromode); 
  fscanf(fp,"%d\n",&mesh.dofluvial); 
  fscanf(fp,"%d\n",&mesh.doice); 
  fscanf(fp,"%d\n",&mesh.dosliding); 
  fscanf(fp,"%d\n",&mesh.doperiglacial); 
  fscanf(fp,"%d\n",&mesh.doglacialhydrology); 
  fscanf(fp,"%d\n",&mesh.doglacialerosion); 
  fscanf(fp,"%d\n",&mesh.dohillslope); 
  fscanf(fp,"%d\n",&mesh.dohillslopeerosion); 
  fscanf(fp,"%d\n",&mesh.doweathering); 
  fscanf(fp,"%d\n",&mesh.dolandslide); 
  fscanf(fp,"%d\n",&mesh.doavalance); 
  fscanf(fp,"%d\n",&mesh.doorographic); 
  fscanf(fp,"%d\n",&mesh.doisostasy); 
  fscanf(fp,"%d\n",&mesh.docelldata); 
  fscanf(fp,"%d\n",&mesh.dosediment); 
  fscanf(fp,"%d\n",&mesh.doglacialsedi); 
  fscanf(fp,"%d\n",&mesh.doparticles); 
  fscanf(fp,"%d\n",&mesh.dodeposit); 
  fscanf(fp,"%d\n",&mesh.dodebrisablation); 
  fclose(fp);

  /*number of cells*/
  mesh.nc = mesh.nx*mesh.ny;

  return(mesh);

}

iproptype read_iprop(char *path,char file[200],FILE *fp)
{

  iproptype iprop;

  /*read from file mesh.input*/
  sprintf(file,"%s/input/iprop.input",path); /*filename*/
  if ((fp = fopen(file,"rt")) == NULL) {
    printf("Cannot open file : %s/input/iprop.input\n",path);
    exit(1);
  }
  fscanf(fp,"%lf\n",&iprop.gamma); 
  fscanf(fp,"%lf\n",&iprop.gamma0); 
  fscanf(fp,"%lf\n",&iprop.Cs); 
  fscanf(fp,"%lf\n",&iprop.latentheat); 
  fscanf(fp,"%lf\n",&iprop.ki); 
  fscanf(fp,"%lf\n",&iprop.rho); 
  fscanf(fp,"%lf\n",&iprop.cp); 
  fscanf(fp,"%lf\n",&iprop.ifac); 
  fscanf(fp,"%lf\n",&iprop.sfac); 
  fscanf(fp,"%lf\n",&iprop.vbfac); 
  fscanf(fp,"%d\n",&iprop.maxitt_v); 
  fscanf(fp,"%d\n",&iprop.maxitt_s); 
  fscanf(fp,"%lf\n",&iprop.C); 
  fscanf(fp,"%lf\n",&iprop.L0); 
  fscanf(fp,"%lf\n",&iprop.minbeta); 
  fscanf(fp,"%lf\n",&iprop.sedifac); 
  fscanf(fp,"%lf\n",&iprop.Kq); 
  fscanf(fp,"%lf\n",&iprop.Ka); 
  fscanf(fp,"%lf\n",&iprop.ap); 
  fscanf(fp,"%lf\n",&iprop.minefac); 
  fscanf(fp,"%lf\n",&iprop.ksg); 
  fclose(fp);

  return(iprop);

}

mproptype read_mprop(char *path,char file[200],FILE *fp)
{

  int i,j;
  mproptype mprop;


  /*read from file mprop.input*/
  sprintf(file,"%s/input/mprop.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/mprop.input\n",path);
    exit(1);
  }
  fread(&mprop.mtype,sizeof(int),1,fp);
  fread(&mprop.avaslope,sizeof(double),1,fp);
  fread(&mprop.avacurv,sizeof(double),1,fp);
  fread(&mprop.maxacc,sizeof(double),1,fp);
  fread(&mprop.dhice,sizeof(double),1,fp);
  fread(&mprop.lrate,sizeof(double),1,fp);
  fread(&mprop.qb,sizeof(double),1,fp);
  fread(&mprop.Ldebris,sizeof(double),1,fp);
  fread(&mprop.nTemp,sizeof(int),1,fp);
  mprop.Temp = malloc(2*sizeof(double*)); for (i=0;i<2;i++) mprop.Temp[i] = malloc(mprop.nTemp*sizeof(double));
  for (j=0;j<mprop.nTemp;j++) for (i=0;i<2;i++) fread(&mprop.Temp[i][j],sizeof(double),1,fp);
  fread(&mprop.nMrate_h,sizeof(int),1,fp);
  mprop.Mrate_h = malloc(2*sizeof(double*)); for (i=0;i<2;i++) mprop.Mrate_h[i] = malloc(mprop.nMrate_h*sizeof(double));
  for (j=0;j<mprop.nMrate_h;j++) for (i=0;i<2;i++) fread(&mprop.Mrate_h[i][j],sizeof(double),1,fp);
  fread(&mprop.nMrate_T,sizeof(int),1,fp);
  mprop.Mrate_T = malloc(3*sizeof(double*)); for (i=0;i<3;i++) mprop.Mrate_T[i] = malloc(mprop.nMrate_T*sizeof(double));
  for (j=0;j<mprop.nMrate_T;j++) for (i=0;i<3;i++) fread(&mprop.Mrate_T[i][j],sizeof(double),1,fp);
  fread(&mprop.dL,sizeof(double),1,fp);
  fread(&mprop.precip,sizeof(double),1,fp);
  fread(&mprop.cp,sizeof(double),1,fp);
  fread(&mprop.ch,sizeof(double),1,fp);
  fclose(fp);

  return(mprop);

}

hwproptype read_hwprop(char *path,char file[200],FILE *fp)
{

  hwproptype hwprop;


  /*read from file hwprop.input*/
  sprintf(file,"%s/input/hwprop.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/hwprop.input\n",path);
    exit(1);
  }
  fread(&hwprop.a2w,sizeof(double),1,fp);
  fread(&hwprop.tscale,sizeof(double),1,fp);
  fread(&hwprop.S0,sizeof(double),1,fp);
  fread(&hwprop.A0,sizeof(double),1,fp);
  fread(&hwprop.kh,sizeof(double),1,fp);
  fread(&hwprop.h0,sizeof(double),1,fp);
  fread(&hwprop.ds,sizeof(double),1,fp);
  fread(&hwprop.ss0,sizeof(double),1,fp);
  fread(&hwprop.kmin,sizeof(double),1,fp);
  fread(&hwprop.alpha,sizeof(double),1,fp);  
  fread(&hwprop.Ls,sizeof(double),1,fp);
  fread(&hwprop.Lc,sizeof(double),1,fp);
  fread(&hwprop.minqw,sizeof(double),1,fp);  
  fread(&hwprop.dtw,sizeof(double),1,fp);
  fclose(fp);

  return(hwprop);

}

hproptype read_hprop(char *path,char file[200],FILE *fp)
{

  hproptype hprop;


  /*read from file hprop.input*/
  sprintf(file,"%s/input/hprop.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/hprop.input\n",path);
    exit(1);
  }
  fread(&hprop.Ks,sizeof(double),1,fp);
  fread(&hprop.sc,sizeof(double),1,fp);
  fread(&hprop.Ke,sizeof(double),1,fp);
  fread(&hprop.Kw,sizeof(double),1,fp);
  fread(&hprop.Ls,sizeof(double),1,fp);
  fread(&hprop.gamma,sizeof(double),1,fp);
  fread(&hprop.Nc,sizeof(int),1,fp);
  fread(&hprop.maxsedi,sizeof(double),1,fp);
  fclose(fp);

  return(hprop);

}

fproptype read_fprop(char *path,char file[200],FILE *fp)
{

  fproptype fprop;


  /*read from file hprop.input*/
  sprintf(file,"%s/input/fprop.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/fprop.input\n",path);
    exit(1);
  }
  fread(&fprop.pr,sizeof(double),1,fp);
  fread(&fprop.rho_s,sizeof(double),1,fp);
  fread(&fprop.Dg,sizeof(double),1,fp);
  fread(&fprop.kw,sizeof(double),1,fp);
  fread(&fprop.tau_c,sizeof(double),1,fp);
  fread(&fprop.Kt,sizeof(double),1,fp);
  fread(&fprop.Ke,sizeof(double),1,fp);

  fclose(fp);

  return(fprop);

}

flexproptype read_flexprop(char *path,char file[200],FILE *fp)
{

  flexproptype flexprop;


  /*read from file flexprop.input*/
  sprintf(file,"%s/input/flexprop.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/flexprop.input\n",path);
    exit(1);
  }
  fread(&flexprop.Te,sizeof(double),1,fp);
  fread(&flexprop.rho_r,sizeof(double),1,fp);
  fread(&flexprop.rho_s,sizeof(double),1,fp);
  fread(&flexprop.rho_a,sizeof(double),1,fp);
  fread(&flexprop.ncall,sizeof(int),1,fp);

  fclose(fp);

  return(flexprop);

}
 

parproptype read_parprop(char *path,char file[200],FILE *fp)
{

  parproptype parprop;


  /*read from file flexprop.input*/
  sprintf(file,"%s/input/parprop.input",path); /*filename*/
  if ((fp = fopen(file,"rt")) == NULL) {
    printf("Cannot open file : %s/input/parprop.input\n",path);
    exit(1);
  } 
  fscanf(fp,"%ld\n",&parprop.npmax); 
  fscanf(fp,"%lf\n",&parprop.minice); 
  fscanf(fp,"%lf\n",&parprop.minsedi); 
  fscanf(fp,"%lf\n",&parprop.maxsedi); 
  fscanf(fp,"%d\n",&parprop.maxpm); 
  fscanf(fp,"%d\n",&parprop.maxp); 
  fscanf(fp,"%d\n",&parprop.minpm); 
  fscanf(fp,"%lf\n",&parprop.minsedim); 
  fclose(fp);

  return(parprop);

}

pgproptype read_pgprop(char *path,char file[200],FILE *fp)
{

  int i,j;

  pgproptype pgprop;

  /*read periglacial properties*/
  sprintf(file,"%s/input/periglacial.input",path); /*filename*/
  if ((fp = fopen(file,"rb")) == NULL) {
    printf("Cannot open file : %s/input/periglacial.input\n",path);
    exit(1);
  }
  fread(&pgprop.nHs,sizeof(int),1,fp);
  fread(&pgprop.nT0,sizeof(int),1,fp);

  /*allocate*/
  pgprop.Hsv = malloc(pgprop.nHs*sizeof(double)); 
  pgprop.T0v = malloc(pgprop.nT0*sizeof(double)); 
  pgprop.CiT = malloc(pgprop.nT0*sizeof(double)); 
  pgprop.Ci = malloc(pgprop.nHs*sizeof(double*)); for (i=0;i<pgprop.nHs;i++) pgprop.Ci[i] = malloc(pgprop.nT0*sizeof(double));
  pgprop.Tr = malloc(pgprop.nHs*sizeof(double*)); for (i=0;i<pgprop.nHs;i++) pgprop.Tr[i] = malloc(pgprop.nT0*sizeof(double));

  for (i=0;i<pgprop.nHs;i++) fread(&pgprop.Hsv[i],sizeof(double),1,fp);
  for (j=0;j<pgprop.nT0;j++) fread(&pgprop.T0v[j],sizeof(double),1,fp);
  for (j=0;j<pgprop.nT0;j++) for (i=0;i<pgprop.nHs;i++) fread(&pgprop.Ci[i][j],sizeof(double),1,fp);
  for (j=0;j<pgprop.nT0;j++) for (i=0;i<pgprop.nHs;i++) fread(&pgprop.Tr[i][j],sizeof(double),1,fp);
  for (j=0;j<pgprop.nT0;j++) fread(&pgprop.CiT[j],sizeof(double),1,fp);
  fread(&pgprop.rho_b,sizeof(double),1,fp);
  fread(&pgprop.rho_s,sizeof(double),1,fp);
  fread(&pgprop.Ke,sizeof(double),1,fp);
  fread(&pgprop.Kt,sizeof(double),1,fp);
  fread(&pgprop.maxsedi,sizeof(double),1,fp);
  fread(&pgprop.maxice,sizeof(double),1,fp);
  fread(&pgprop.minslope,sizeof(double),1,fp);
  fread(&pgprop.minsedi,sizeof(double),1,fp);
  fclose(fp);

  return(pgprop);

}



void write_output(celltype **cells,hptype **hp,vptype **vp,meshtype mesh,char *path,char file[200],FILE *fp,long step)
{

  int i,j,k;
  double vx,vy;
  int nx = mesh.nx;
  int ny = mesh.ny;

  /*create output file*/
  sprintf(file,"%s/output/output%ld.dat",path,step);
  if ((fp = fopen(file,"wb")) == NULL) {
    printf("could not open file for writing\n");
    exit(1);
  }
  else {
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].x,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].y,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].bed,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].ice,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].bslope,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].hydro,sizeof(int),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].slidingslope,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].ssedi,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].bsedi,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].msedi,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].sN10,sizeof(double),1,fp); 
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].np,sizeof(int),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].sedi,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].tn,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].te,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].ts,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].sliding,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].deformation,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].td2,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Mb,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].quarrying,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].accrate,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Ms,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].oro,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].precip,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].weathering,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].abrasion,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].isostasy,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].periglacial_erosion,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Ts,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].hillslope_erosion,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Ta,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Tb,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].sfac,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].margin,sizeof(int),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Ac,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Kgw,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].qw,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].afac,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].water,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Pw,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].SLf,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].mw,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].sxx,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].sxy,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].syy,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].dsxxdx,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].dsxydx,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].dsyydx,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].dsxxdy,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].dsxydy,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].dsyydy,sizeof(double),1,fp);
    for (j=0;j<nx+1;j++) for (i=0;i<ny;i++) fwrite(&hp[i][j].vx,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny+1;i++) fwrite(&vp[i][j].vy,sizeof(double),1,fp);
    for (j=0;j<nx+1;j++) for (i=0;i<ny;i++) fwrite(&hp[i][j].vresx,sizeof(double),1,fp);
    for (j=0;j<nx;j++) for (i=0;i<ny+1;i++) fwrite(&vp[i][j].vresy,sizeof(double),1,fp);
    /*for (k=0;k<20;k++) for (j=0;j<nx;j++) for (i=0;i<ny;i++) fwrite(&cells[i][j].Vs[k],sizeof(double),1,fp);*/
    fclose(fp);
  }
    
}/*write_output*/
 

void write_particles(particletype *pp,long npa,long *pactive,char *path,char file[200],FILE *fp,long step)
{

  int i;

  /*create output file*/
  sprintf(file,"%s/output/particles%ld.dat",path,step);
  if ((fp = fopen(file,"wb")) == NULL) {
    printf("could not open file for writing\n");
    exit(1);
  }
  else {
    fwrite(&npa,sizeof(int),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].x,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].y,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].z,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].bf,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].sedi,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].N10,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].bx,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].by,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].birthday,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].age,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].dl,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].vx,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].vy,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].vz,sizeof(double),1,fp);
    for (i=0;i<npa;i++) fwrite(&pp[pactive[i]].erate,sizeof(double),1,fp);

    fclose(fp);
  }
    
}/*write_output*/

